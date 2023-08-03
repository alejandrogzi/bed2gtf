use regex::Regex;

use std::string::String;
use std::path::PathBuf;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::collections::{HashMap, HashSet};
use std::io::Write;

mod dp;
use dp::Dependencies;

mod sh;
use sh::Shell;



/// bed2gtf currently uses a combined approach between a Rust implementation
/// and UCSC binaries. -> this is a temporary solution until a better Rust
/// implementation is developed
/// run_binary automates the download of binaries through the dp module and 
/// runs them using the sh module.
/// cmd -> bedToGenePred {args} stdout | genePredToGtf stdin {output}
pub fn run_binary(bed: PathBuf) -> PathBuf {
    let gtf = bed.parent().unwrap().join("temp.gtf");
    let binaries = Dependencies::get();
    let bed_to_gene_pred = &binaries[0];
    let gene_pred_to_gtf = &binaries[1];

    let to_gene_pred: String = format!("{} {} stdout", bed_to_gene_pred.to_string_lossy(), bed.display());
    let to_gtf: String = format!("{} file stdin {}", gene_pred_to_gtf.to_string_lossy(), gtf.to_string_lossy());
    let cmd: String = format!("{} | {}", to_gene_pred, to_gtf);

    let _ = Shell::run(&cmd);
    return gtf;
}



/// get_isoforms stores transcript:gene pairs in a HashMap
/// to be used by the mapper function.
pub fn get_isoforms(path: PathBuf) -> Result<HashMap<String, String>, io::Error> {
    let file: File = File::open(path).unwrap();
    let reader: BufReader<File> = BufReader::new(file);
    let mut isoforms: HashMap<String, String> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let content: Vec<&str> = line
        .split("\t")
        .collect();

        let gene: &str = content[0];
        let transcript: &str = content[1];
        isoforms.insert(transcript.to_string(), gene.to_string());
    }

    return Ok(isoforms);
} 



/// fix_gtf replaces the gene_id in the GTF file with the gene_id 
/// in the isoforms mapping. 
/// this also could handle the mapping of gene_biotype features but
/// will be experimental for some time -> need to test efficiency cost.
pub fn fix_gtf(gtf: PathBuf, isoforms: PathBuf) -> Result<PathBuf, io::Error> {
    let output_file_path = gtf.with_extension("fixed");
    let isoforms_mapping: HashMap<String, String> = get_isoforms(isoforms)?;
    let gtf_file = File::open(gtf)?;
    let reader = BufReader::new(gtf_file);
    let mut new_gtf_file = File::create(&output_file_path)?;
    let re: Regex = Regex::new(r#"gene_id\s+"([^"]+)""#).unwrap();

    for line in reader.lines() {
        let line = line?; // Propagate any read errors
        let mut fields: Vec<&str> = line.split('\t').collect();

        // Check if there are at least 9 fields in the line
        if fields.len() < 9 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid GTF format: Less than 9 columns in the line.",
            ));
        }

        let mut attributes: Vec<&str> = fields[8].split(';').collect(); // ["gene_id \"ENSG00000223972.5\"", ...]
        let mut gene: Vec<&str> = attributes[0].split('"').collect(); // ["gene_id ", "ENSG00000223972.5", ""]

        // Check if the line contains the gene_id information
        if let Some(captures) = re.captures(&line) {
            let transcript = captures
                .get(1)
                .map_or("", |m| m.as_str())
                .split('.')
                .next()
                .unwrap_or("");

            // Look up the gene_id for the current transcript in the isoforms mapping
            if let Some(result) = isoforms_mapping.get(transcript) {
                gene[1] = result;
                let fixed_gene: String = gene.join("\"");

                attributes[0] = fixed_gene.as_str();
                let fixed_attributes: String = attributes.join(";");

                fields[8] = fixed_attributes.as_str();
                let fixed_metadata: String = fields.join("\t");

                // Write the fixed line to the output file
                writeln!(new_gtf_file, "{}", fixed_metadata)?;

            } else {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("Gene not found for transcript: {}", transcript),
                ));
            }
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("gene_id not found in line: {}", line),
            ));
        }
    }
    Ok(output_file_path)
}



/// inserts a 'gene' feature in GTF file by grabbing the 
/// first transcript-featured appearance of a gene_id.
pub fn insert_gene(gtf: PathBuf) -> Result<PathBuf, io::Error> {
    let output_file_path: PathBuf = gtf.parent().unwrap().join("output.gtf");
    let mut new_gtf_file: File = File::create(&output_file_path)?;
    let gtf_file: File = File::open(gtf)?;
    let reader: BufReader<File> = BufReader::new(gtf_file);
    let mut temp: HashSet<String> = HashSet::new();


    for line in reader.lines() {
        if let Ok(line) = line {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            let feature_type = fields[2];

            if feature_type != "transcript" {
                writeln!(new_gtf_file, "{}", line)?;
            } else {
                let gene_id: &str = fields[8]
                    .split(';').find(|s| s.starts_with("gene_id"))
                    .unwrap();
                let gene: &str = gene_id.split('"').nth(1).unwrap();
                
                if !temp.contains(gene) {
                    temp.insert(gene.to_string());

                    let mut base = fields[0..9].to_owned();
                    let new_attribute = &format!("{}; gene_biotype \"protein_coding\";",gene_id);
                    base[2] = "gene";
                    base[8] = new_attribute;

                    let newline = base.join("\t");
                    writeln!(new_gtf_file, "{}", newline)?;
                }
                writeln!(new_gtf_file, "{}", line)?;
            }
        }
    }
    Ok(output_file_path)
}



pub fn clean_up(path: PathBuf) -> Result<(), io::Error> {
    let mut in_path = path.clone();
    in_path.pop();

    let stem = path.file_stem().unwrap().to_str().unwrap();

    let files = std::fs::read_dir(in_path)?;
    for file in files {
        let file = file.unwrap().path();
        if file.is_file() {
            if file.file_stem().unwrap() != stem {
                continue
            } else {
                std::fs::remove_file(file)?;
            }
        }
    }
    Ok(())
}


pub fn bed2gtf(bed: PathBuf, isoforms: PathBuf) -> PathBuf {
    let input: PathBuf = run_binary(bed);
    let gtf: PathBuf = {
        let gtf: PathBuf = fix_gtf(input.clone().into(), isoforms).unwrap();
        insert_gene(gtf).unwrap()
    };
    clean_up(input).unwrap();
    gtf
}