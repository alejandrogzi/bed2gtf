use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::fs::File;
use std::cmp::{max, min};
use std::time::Instant;
use std::error::Error;

use colored::Colorize;
use peak_alloc::PeakAlloc;

use log::Level;

use indoc::indoc;

mod bed;
use bed::BedRecord;

mod codon;
use codon::Codon;


const SOURCE: &str = "bed2gtf";

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;


/// Store the isoforms in a HashMap in the form: isoform -> gene.
/// Provides a fast access to the gene name given an isoform name.
fn get_isoforms(path: PathBuf) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let file: File = File::open(path).unwrap();
    let reader: BufReader<File> = BufReader::new(file);
    let mut isoforms: HashMap<String, String> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let content: Vec<&str> = line
        .split("\t")
        .collect();

        let gene: &str = content[0];
        let isoform: &str = content[1];
        isoforms.insert(isoform.to_string(), gene.to_string());
    }

    return Ok(isoforms);
} 



/// Get the coordinates of the first codon (start/stop).
/// If not in frame, return an empty codon.
fn find_first_codon(record: &BedRecord) -> Codon {
    let mut codon = Codon::new();
    let mut exon = 0;

    for k in 0..record.get_exon_frames().len() {
        if record.get_exon_frames()[k] >= 0 {
            exon = k;
            break;
        } else {
            return codon;
        }
    }

    let cds_exon_start = max(record.exon_start()[exon], record.cds_start());
    let cds_exon_end = min(record.exon_end()[exon], record.cds_end());

    let frame = if record.strand() == "+" {
        record.get_exon_frames()[exon]
    } else {
        (record.get_exon_frames()[exon] + (cds_exon_end - cds_exon_start)) % 3
    };

    if frame != 0 {
        return codon;
    }

    codon.start = record.cds_start();
    codon.end = record.cds_start() + (record.cds_end() - record.cds_start()).min(3);
    codon.index = exon as i32;

    if codon.end - codon.start < 3 {
        exon = exon + 1;
        if exon == record.exon_count() as usize {
            return codon;
        };
        let need = 3 - (codon.end - codon.start);
        if (record.cds_end() - record.cds_start()) < need {
            return codon;
        }
        codon.start2 = record.cds_start();
        codon.end2 = record.cds_start() + need;
    }
    codon
}



/// Get the coordinates of the last codon (start/stop).
/// If not in frame, return an empty codon.
fn find_last_codon(record: &BedRecord) -> Codon {
    let mut codon = Codon::new();
    let mut exon = 0;

    for k in (0..record.get_exon_frames().len()).step_by(1).rev() {
        if record.get_exon_frames()[k] >= 0 {
            exon = k;
            break;
        } else {
            return codon;
        }
    }

    let cds_exon_start = max(record.exon_start()[exon], record.cds_start());
    let cds_exon_end = min(record.exon_end()[exon], record.cds_end());

    let frame = if record.strand() == "-" {
        record.get_exon_frames()[exon]
    } else {
        (record.get_exon_frames()[exon] + (cds_exon_end - cds_exon_start)) % 3
    };

    if frame != 0 {
        return codon;
    }

    codon.start = max(record.cds_start(), record.cds_end() - 3);
    codon.end = record.cds_end();
    codon.index = exon as i32;

    if codon.end - codon.start < 3 {
        exon = exon + 1;
        if exon == record.exon_count() as usize {
            return codon;
        };
    
        let need = 3 - (codon.end - codon.start);
        if (record.cds_end() - record.cds_start()) < need {
            return codon;
        }
        codon.start2 = record.cds_start();
        codon.end2 = record.cds_start() + need;
    }
    codon
}



/// Check if all the bases of a codon are defined.
fn codon_complete(codon: &Codon) -> bool {
    ((codon.end - codon.start) + (codon.end2 - codon.start2)) == 3
}



/// Check if a given coordinate is within exon boundaries.
fn in_exon(record: &BedRecord, pos:i32, exon: usize) -> bool {
    (record.exon_start()[exon] <= pos) && (pos <= record.exon_end()[exon])
}



/// Move a position in an exon by a given distance, which is positive
/// to move forward and negative to move backwards. Introns are not
/// considered. If can't move distance and stay within exon boundaries,
/// panic.
fn move_pos(record: &BedRecord, pos: i32, dist: i32) -> i32 {
    let mut pos = pos;
    assert!(record.tx_start() <= pos && pos <= record.tx_end());
    let mut exon: Option<i16> = None;
    for i in 0..record.exon_count() {
        if in_exon(record, pos, i as usize) {
            exon = Some(i);
            break;
        }
    } 

    if exon.is_none() {
        panic!("Position {} not in exons", pos);
    }

    let mut steps = dist.abs();
    let direction = if dist >= 0 { 1 } else { -1 };

    while (0..record.exon_count()).contains(&(exon.unwrap())) && steps > 0 {
        if in_exon(record, pos + direction, exon.unwrap() as usize) {
            pos += direction;
            steps -= 1;
        } else if direction >= 0 {
            exon = Some(exon.unwrap() + 1);
            if let Some(ex) = exon {
                if (ex as usize) < record.exon_count() as usize {
                    pos = record.exon_start()[ex as usize];
                }
            }
        } else {
            exon = Some(exon.unwrap() - 1);
            if let Some(ex) = exon {
                if ex >= 0 {
                    pos = record.exon_end()[ex as usize] - 1;
                    steps -= 1;
                }
            }
        }
    }

    if steps > 0 {
        panic!("can't move {} by {}", pos, dist);
    }

    pos
} 



/// Build a "gene" feature line for a given group of transcripts.
/// Each line is unique for a given group.
fn build_gene_line(gene_name: &str, record: &BedRecord, file: &mut File) {
    assert!(gene_name.len() > 0);
    let gene_line = format!("{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tgene_id \"{}\";\n",
        record.chrom(),
        SOURCE,
        record.tx_start() + 1,
        record.tx_end(),
        record.strand(),
        gene_name
    );
    file.write_all(gene_line.as_bytes()).unwrap();
}



/// Build a GTF line for a given feature.
fn build_gtf_line(record: &BedRecord, gene_name: &str, gene_type: &str, exon_start: i32, exon_end: i32, frame: i32, exon: i16, file: &mut File) {
    assert!(record.tx_start() < record.tx_end());

    let phase = if frame < 0 {
        "."
    } else if frame == 0 {
        "0"
    } else if frame == 1 {
        "2"
    } else {
        "1"
    };

    let mut gtf_line = format!("{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t",
        record.chrom(),
        SOURCE,
        gene_type,
        exon_start + 1,
        exon_end,
        record.strand(),
        phase,
    );

    gtf_line += &format!("gene_id \"{}\"; ", gene_name);
    gtf_line += &format!("transcript_id \"{}\";", record.name());
    if exon >= 0 {
        if record.strand() == "-" {
            gtf_line += &format!(" exon_number \"{}\";", record.exon_count() - exon);
            gtf_line += &format!(" exon_id \"{}.{}\";", record.name(), record.exon_count() - exon);
        } else {
            gtf_line += &format!(" exon_number \"{}\";", exon + 1);
            gtf_line += &format!(" exon_id \"{}.{}\";", record.name(), exon + 1);
        }
    }
    gtf_line += "\n";
    let _ = file.write_all(gtf_line.as_bytes());
}



/// Write the features of a given exon, including UTRs and CDS.
fn write_features(i: usize, record: &BedRecord, gene_name: &str, first_utr_end: i32, cds_start: i32, cds_end: i32, last_utr_start: i32, frame: i32, file: &mut File) {
    let exon_start = record.exon_start()[i];
    let exon_end = record.exon_end()[i];

    if exon_start < first_utr_end {
        let end = min(exon_end, first_utr_end);
        let utr_type = if record.strand() == "+" { "five_prime_utr" } else { "three_prime_utr" };
        build_gtf_line(record, gene_name, utr_type, exon_start, end, frame, -1, file);
    }

    if record.cds_start() < exon_end && exon_start < record.cds_end() {
        let start = max(exon_start, cds_start);
        let end = min(exon_end, cds_end);
        build_gtf_line(record, gene_name, "CDS", start, end, frame, i as i16, file);
    }

    if exon_end > last_utr_start {
        let start = max(exon_start, last_utr_start);
        let utr_type = if record.strand() == "+" { "three_prime_utr" } else { "five_prime_utr" };
        build_gtf_line(record, gene_name, utr_type, start, exon_end, frame, -1, file);
    }
}



/// Write the codon features (start/stop) for a given exon.
fn write_codon(record: &BedRecord, gene_name: &str, gene_type: &str, codon: Codon, file: &mut File) {
    build_gtf_line(record, gene_name, gene_type, codon.start, codon.end, 0, codon.index as i16, file);

    if codon.start2 < codon.end2 {
        build_gtf_line(record, gene_name, gene_type, codon.start, codon.end, codon.start2, (codon.end - codon.start) as i16, file);
    }
}



/// Convert a BED record to a GTF record.
fn to_gtf(record: &BedRecord, isoforms: &HashMap<String, String>, file: &mut File, gene_line: bool) {
    let gene_name = isoforms.get(record.name()).unwrap();
    let first_codon = find_first_codon(record);
    let last_codon = find_last_codon(record);

    let first_utr_end = record.cds_start();
    let last_utr_start = record.cds_end();

    let cds_end: i32 = if record.strand() == "+" && codon_complete(&last_codon) {
        move_pos(record, last_codon.end, -3)
    } else {
        record.cds_end()
    };

    let cds_start = if record.strand() == "-" && codon_complete(&first_codon) {
        move_pos(record, first_codon.start, 3)
    } else {
        record.cds_start()
    };

    if gene_line {
        let _ = build_gene_line(gene_name, record, file);
    };

    let _ = build_gtf_line(record, gene_name, "transcript", record.tx_start(), record.tx_end(), -1, -1, file);

    for i in 0..record.exon_count() as usize {
        build_gtf_line(record, gene_name, "exon", record.exon_start()[i], record.exon_end()[i], -1, i as i16, file);
        if cds_start < cds_end {
            write_features(i, record, gene_name, first_utr_end, cds_start, cds_end, last_utr_start, record.get_exon_frames()[i], file);
        }
    }

    if record.strand() != "-" {
        if codon_complete(&first_codon) {
            write_codon(record, gene_name, "start_codon", first_codon, file);
        }
        if codon_complete(&last_codon) {
            write_codon(record, gene_name, "stop_codon", last_codon, file);
        }
    } else {
        if codon_complete(&last_codon) {
            write_codon(record, gene_name, "start_codon", last_codon, file);
        }
        if codon_complete(&first_codon) {
            write_codon(record, gene_name, "stop_codon", first_codon, file);
        }
    }
}



/// Convert a BED file to a GTF file.
/// ```
/// use bed2gtf::bed2gtf;
/// bed2gtf("input.bed", "isoforms.txt", "output.gtf");
/// ```
pub fn bed2gtf(input: &String, isoforms: &String, output: &String) -> Result<(), Box<dyn Error>> {

    msg();
    simple_logger::init_with_level(Level::Info)?;

    let start = Instant::now();
    let bedfile = File::open(PathBuf::from(input)).unwrap();
    let reader = BufReader::new(bedfile);
    
    let isoforms = get_isoforms(isoforms.into()).unwrap();
    let mut output = File::create(PathBuf::from(output)).unwrap();
    let mut seen_genes: HashSet<String> = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let record = BedRecord::new(fields);

        let key = isoforms.get(record.name()).unwrap();
        if !seen_genes.contains(key) {
            seen_genes.insert(key.to_string());
            let _ = to_gtf(&record, &isoforms, &mut output, true);
        } else {
            let _ = to_gtf(&record, &isoforms, &mut output, false);
        }
    }

    let peak_mem = PEAK_ALLOC.peak_usage_as_mb();

    log::info!("Memory usage: {} MB", peak_mem);
    log::info!("Elapsed: {:.4?}", start.elapsed().as_secs_f32());

    Ok(())
}



fn msg() {
    println!("{}\n{}",
        "\n##### BED2GTF #####".bright_cyan().bold(),
        indoc!("A fast BED-to-GTF converter written in Rust.
        Repository: https://github.com/alejandrogzi/bed2gtf
        Feel free to contact the developer if any issue/suggest/bug is found.
        "));
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn make_bed_test_file(content: &str) -> Result<String, Box<dyn Error>> {
        let mut file = File::create("test.bed")?;
        file.write_all(content.as_bytes())?;

        Ok("test.bed".to_string())
    }

    fn make_isoforms_test_file(content: &str) -> Result<String, Box<dyn Error>> {
        let mut file = File::create("test.isoforms")?;
        file.write_all(content.as_bytes())?;

        Ok("test.isoforms".to_string())
    }

    #[test]
    fn main_test() {
        let input = make_bed_test_file("chr15\t81000922\t81005788\tENST00000267984\t0\t+\t81002271\t81003360\t0\t1\t4866,\t0,").unwrap();

        let isoforms = make_isoforms_test_file(indoc!("ENSG00000140406\tENST00000267984
        ENSG00000140545\tENST00000560937
        ENSG00000082438\tENST00000495084")).unwrap();

        let output = "out.gtf".to_string();

        let _ = bed2gtf(&input, &isoforms, &output).unwrap();
        let mut file = File::open(output.clone()).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();

        let converted_content = indoc!("chr15\tbed2gtf\tgene\t81000923\t81005788\t.\t+\t.\tgene_id \"ENSG00000140406\";
        chr15\tbed2gtf\ttranscript\t81000923\t81005788\t.\t+\t.\tgene_id \"ENSG00000140406\"; transcript_id \"ENST00000267984\";
        chr15\tbed2gtf\texon\t81000923	81005788\t.\t+\t.\tgene_id \"ENSG00000140406\"; transcript_id \"ENST00000267984\"; exon_number \"1\"; exon_id \"ENST00000267984.1\";
        chr15\tbed2gtf\tfive_prime_utr\t81000923\t81002271\t.\t+\t0\tgene_id \"ENSG00000140406\"; transcript_id \"ENST00000267984\";
        chr15\tbed2gtf\tCDS	81002272\t81003357\t.\t+\t0\tgene_id \"ENSG00000140406\"; transcript_id \"ENST00000267984\"; exon_number \"1\"; exon_id \"ENST00000267984.1\";
        chr15\tbed2gtf\tthree_prime_utr\t81003361\t81005788\t.\t+\t0\tgene_id \"ENSG00000140406\"; transcript_id \"ENST00000267984\";
        chr15\tbed2gtf\tstart_codon\t81002272\t81002274\t.\t+\t0\tgene_id \"ENSG00000140406\"; transcript_id \"ENST00000267984\"; exon_number \"1\"; exon_id \"ENST00000267984.1\";
        chr15\tbed2gtf\tstop_codon\t81003358\t81003360\t.\t+\t0\tgene_id \"ENSG00000140406\"; transcript_id \"ENST00000267984\"; exon_number \"1\"; exon_id \"ENST00000267984.1\";\n");

        assert_eq!(contents, converted_content);

        teardown(vec![input, isoforms, output]);
    }

    fn teardown(files: Vec<String>) {
        for file in files {
            std::fs::remove_file(file).unwrap();
        }
    }
}