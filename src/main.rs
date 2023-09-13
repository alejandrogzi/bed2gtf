use bed2gtf::bed2gtf;

use clap::{Arg, Command, ArgMatches};

use colored::Colorize;

use std::string::String;
use std::error::Error;



fn main() {
    let matches = Command::new("bed2gtf")
        .version("1.5.0")
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
        .arg(Arg::new("output")
            .index(3)
            .required(true)
            .value_name("OUTPUT")
            .help("Output file name"))
        .get_matches();

    if let Some(err) = run(matches).err() {
        eprintln!("{} {}", 
                "Error:".bright_red().bold(),
                err);
        std::process::exit(1);
    }
}


fn run(matches: ArgMatches) -> Result<(), Box<dyn Error>> {
    let bed: &String = matches.get_one("bed").unwrap();
    let isoforms: &String = matches.get_one("isoforms").unwrap();
    let output: &String = matches.get_one("output").unwrap();

    let _ = bed2gtf(bed, isoforms, output);

    println!("{} {}", 
    "Success:".bright_green().bold(),
    "BED file converted successfully!");

    Ok(())
}