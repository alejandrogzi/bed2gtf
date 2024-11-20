use crate::bed::BedRecord;

use chrono::Datelike;
use colored::Colorize;
use flate2::read::GzDecoder;
use indoc::indoc;
use rayon::prelude::*;

use std::collections::HashMap;
use std::error::Error;
use std::fmt::Debug;
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};

const SOURCE: &str = "bed2gtf";
const VERSION: &str = env!("CARGO_PKG_VERSION");
const REPOSITORY: &str = env!("CARGO_PKG_REPOSITORY");

pub fn get_isoforms(file: &String) -> HashMap<String, String> {
    let pairs = parallel_hash_rev(file);
    // let rev_pairs = parallel_hash(&file);

    if pairs.len() == 0 {
        println!(
            "{} {}",
            "Fail:".bright_red().bold(),
            "BED file could not be converted. Please check your isoforms file."
        );
        std::process::exit(1);
    }
    // (pairs, rev_pairs)
    pairs
}

pub fn reader(file: &PathBuf) -> io::Result<String> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

pub fn raw<P: AsRef<Path> + Debug>(f: P) -> Result<String, Box<dyn Error>> {
    let mut file = File::open(f)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

pub fn with_gz<P: AsRef<Path> + Debug>(f: P) -> Result<String, Box<dyn Error>> {
    let file = File::open(f)?;
    let mut decoder = GzDecoder::new(BufReader::new(file));

    let mut contents = String::new();
    decoder.read_to_string(&mut contents)?;

    Ok(contents)
}

pub fn parallel_hash<'a>(s: &'a str) -> HashMap<String, String> {
    s.par_lines()
        .filter_map(|line| {
            let mut words = line.split_whitespace();
            if let Some(fw) = words.next() {
                if let Some(sw) = words.next() {
                    return Some((fw.to_owned(), sw.to_owned()));
                }
            }
            // If the line doesn’t have two words, ignore it and return None.
            None
        })
        .collect()
}

pub fn parallel_hash_rev<'a>(s: &'a str) -> HashMap<String, String> {
    s.par_lines()
        .filter_map(|line| {
            let mut words = line.split_whitespace();
            if let Some(fw) = words.next() {
                if let Some(sw) = words.next() {
                    return Some((sw.to_owned(), fw.to_owned()));
                }
            }
            None
        })
        .collect()
}

pub fn parallel_parse<'a>(s: &'a str) -> Result<Vec<BedRecord>, String> {
    let records = s
        .par_lines()
        // .map(|line| BedRecord::parse(line))
        .filter_map(|line| match std::str::from_utf8(line.as_bytes()) {
            Ok(valid_line) => Some(BedRecord::parse(valid_line)),
            Err(_) => {
                eprintln!("Skipping invalid UTF-8 line: {:?}", line);
                None
            }
        })
        .collect::<Result<Vec<BedRecord>, String>>();

    Ok(records?)
}

pub fn custom_par_parse(
    records: &Vec<BedRecord>,
) -> Result<HashMap<String, (String, u32, u32, String)>, &'static str> {
    let gene_coordinates = records
        .into_par_iter()
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, (String, u32, u32, String)>, record| {
                acc.entry(record.name.clone()).or_insert((
                    record.chrom.clone(),
                    record.tx_start,
                    record.tx_end,
                    record.strand.clone(),
                ));
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut a: HashMap<String, (String, u32, u32, String)>, b| {
                for (key, (chrom, start, end, strand)) in b {
                    a.entry(key).or_insert((chrom, start, end, strand));
                }
                a
            },
        );
    Ok(gene_coordinates)
}

pub fn combine_maps_par(
    isoforms: &HashMap<String, String>,
    gene_track: &HashMap<String, (String, u32, u32, String)>,
) -> Vec<(String, String, u32, u32, String, String, String)> {
    let coords = isoforms
        .par_iter()
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, (String, u32, u32, String)>, (transcript, gene)| {
                if let Some(&(ref chrom, start, end, ref strand)) = gene_track.get(transcript) {
                    let entry = acc.entry(gene.clone()).or_insert((
                        chrom.to_string(),
                        start,
                        end,
                        strand.to_string(),
                    ));
                    entry.1 = entry.1.min(start); // Update min start
                    entry.2 = entry.2.max(end); // Update max end
                }
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut a, b| {
                for (gene, (chrom, start, end, strand)) in b {
                    let entry = a.entry(gene).or_insert((chrom, start, end, strand));
                    entry.1 = entry.1.min(start); // Update min start
                    entry.2 = entry.2.max(end); // Update max end
                }
                a
            },
        );

    let lines = coords
        .par_iter()
        .map(|(gene, (chrom, start, end, strand))| {
            (
                chrom.to_string(),
                "gene".to_string(),
                start + 1,
                *end,
                strand.to_string(),
                ".".to_string(),
                format!("gene_id \"{}\";", gene),
            )
        })
        .collect();
    lines
}

pub fn max_mem_usage_mb() -> f64 {
    let rusage = unsafe {
        let mut rusage = std::mem::MaybeUninit::uninit();
        libc::getrusage(libc::RUSAGE_SELF, rusage.as_mut_ptr());
        rusage.assume_init()
    };
    let maxrss = rusage.ru_maxrss as f64;
    if cfg!(target_os = "macos") {
        maxrss / 1024.0 / 1024.0
    } else {
        maxrss / 1024.0
    }
}

pub fn msg() {
    println!(
        "{}\n{}\n{}\n",
        "\n##### BED2GTF #####".bright_cyan().bold(),
        indoc!(
            "A fast BED-to-GTF converter written in Rust.
        Repository: https://github.com/alejandrogzi/bed2gtf
        Feel free to contact the developer if any issue/bug is found."
        ),
        format!("Version: {}", VERSION)
    );
}

pub fn get_date() -> String {
    let now = chrono::Utc::now();
    let year = now.year();
    let month = now.month();
    let day = now.day();

    format!("{}-{}-{}", year, month, day)
}

pub fn comments(file: &mut Box<dyn Write>) {
    let _ = file.write_all(format!("#provider: {}\n", SOURCE).as_bytes());
    let _ = file.write_all(format!("#version: {}\n", VERSION).as_bytes());
    let _ = file.write_all(format!("#contact: {}\n", REPOSITORY).as_bytes());
    let _ = file.write_all(format!("#date: {}\n", get_date()).as_bytes());
}
