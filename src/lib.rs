use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write, Error};
use std::path::PathBuf;
use std::fs::File;
use std::cmp::{max, min};

mod bed;
use bed::BedRecord;

mod codon;
use codon::Codon;


const SOURCE: &str = "bed2gtf";


fn get_isoforms(path: PathBuf) -> Result<HashMap<String, String>, Error> {
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


fn codon_complete(codon: &Codon) -> bool {
    ((codon.end - codon.start) + (codon.end2 - codon.start2)) == 3
}


fn in_exon(record: &BedRecord, pos:i32, exon: usize) -> bool {
    (record.exon_start()[exon] <= pos) && (pos <= record.exon_end()[exon])
}


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


fn build_gene_line(gene_name: &str, record: &BedRecord, file: &mut File) {
    assert!(gene_name.len() > 0);
    let gene_line = format!("{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tgene_id \"{}\";\n",
        record.chrom(),
        SOURCE,
        record.tx_start() + 1,
        record.tx_end(),
        record.strand(),
        gene_name,
    );
    file.write_all(gene_line.as_bytes()).unwrap();
}


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
    gtf_line += &format!("transcript_id \"{}\"; ", record.name());
    if exon >= 0 {
        if record.strand() == "-" {
            gtf_line += &format!("exon_number \"{}\"; ", record.exon_count() - exon);
            gtf_line += &format!("exon_id \"{}.{}\";", record.name(), record.exon_count() - exon);
        } else {
            gtf_line += &format!("exon_number \"{}\"; ", exon + 1);
            gtf_line += &format!("exon_id \"{}.{}\";", record.name(), exon + 1);
        }
    }
    gtf_line += "\n";
    let _ = file.write_all(gtf_line.as_bytes());
}



fn write_features(i: usize, record: &BedRecord, gene_name: &str, first_utr_end: i32, cds_start: i32, cds_end: i32, last_utr_start: i32, frame: i32, file: &mut File) {
    let exon_start = record.exon_start()[i];
    let exon_end = record.exon_end()[i];

    if exon_start < first_utr_end {
        let end = min(exon_end, first_utr_end);
        let utr_type = if record.strand() == "+" { "5UTR" } else { "3UTR" };
        build_gtf_line(record, gene_name, utr_type, exon_start, end, frame, -1, file);
    }

    if record.cds_start() < exon_end && exon_start < record.cds_end() {
        let start = max(exon_start, cds_start);
        let end = min(exon_end, cds_end);
        build_gtf_line(record, gene_name, "CDS", start, end, frame, i as i16, file);
    }

    if exon_end > last_utr_start {
        let start = max(exon_start, last_utr_start);
        let utr_type = if record.strand() == "+" { "3UTR" } else { "5UTR" };
        build_gtf_line(record, gene_name, utr_type, start, exon_end, frame, -1, file);
    }
}


fn write_codon(record: &BedRecord, gene_name: &str, gene_type: &str, codon: Codon, file: &mut File) {
    build_gtf_line(record, gene_name, gene_type, codon.start, codon.end, 0, codon.index as i16, file);

    if codon.start2 < codon.end2 {
        build_gtf_line(record, gene_name, gene_type, codon.start, codon.end, codon.start2, (codon.end - codon.start) as i16, file);
    }
}


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
        build_gene_line(gene_name, record, file)
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



pub fn bed2gtf(input: &String, isoforms: &String, output: &String) {
    let bedfile = File::open(PathBuf::from(input)).unwrap();
    let reader = BufReader::new(bedfile);
    
    let isoforms = get_isoforms(isoforms.into()).unwrap();
    let mut output = File::create(PathBuf::from(output)).unwrap();
    let mut seen_genes: HashSet<String> = HashSet::new();

    for line in reader.lines() {
        let line = line.unwrap();
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
}