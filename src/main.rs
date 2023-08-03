use bed2gtf::*;

use clap::{Arg, Command, ArgMatches};

use colored::Colorize;

use indoc::indoc;

use std::string::String;
use std::path::PathBuf;
use std::error::Error;
use std::time::Instant;



fn main() {
    let matches = Command::new("bed2gtf")
        .version("1.0")
        .author("Alejandro Gonzales-Irribarren <jose.gonzalesdezavala1@unmsm.edu.pe>")
        .about("A fast and memory efficient BED to GTF converter")
        .arg(Arg::new("bed")
            .index(1)
            .required(true)
            .value_name("BED")
            .help("BED file to convert"))
        .arg(Arg::new("isoforms")
            .index(2)
            .required(true)
            .value_name("ISOFORMS")
            .help("Isoforms mapping file"))
        .arg(Arg::new("verbose")
            .index(3)
            .default_value("true")
            .help("Prints verbose output"))
        .get_matches();

    if let Some(err) = run(matches).err() {
        eprintln!("{} {}", 
                "Error:".bright_red().bold(),
                err.to_string().bright_red().bold());
        std::process::exit(1);
    }
}


//TO DO: add a flag to keep the temp files
//TO DO: improve print statements

/// bed2gtf runner; automates whole process
fn run(matches: ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed: &String = matches.get_one("bed").unwrap();
    let isoforms: &String = matches.get_one("isoforms").unwrap();

    let now = Instant::now();
    
    println!("{}", indoc!(
    "\n
    ##### BED2GTF in RUST #####
    A fast and memory efficient BED to GTF converter.\n
    ##### STEP 1: DOWNLOAD DEPENDENCIES #####
    "));

    let gtf: PathBuf = run_binary(bed.into());

    println!("\n##### STEP 2: CONVERTING BED TO GTF #####");
    println!("Fixing gene_ids in the GTF file...");
    let gtf: PathBuf = fix_gtf(gtf.into(), isoforms.into())?;
    println!("Done!\n");

    println!("##### STEP 3: INSERT GENE FEATURES #####");
    println!("Inserting gene features in the GTF file...");
    let out: PathBuf = insert_gene(gtf.clone()).unwrap();
    println!("Done!\n");

    println!("Cleaning up...");
    let _ = clean_up(gtf);

    if matches.contains_id("verbose") {
        println!("results at: {}", out.display());
        let elapsed = now.elapsed();
        println!("time taken: {:.2?}", elapsed);
    }

    Ok(())
}
