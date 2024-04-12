use clap::{self, ArgAction, Parser};
use num_cpus;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Parser, Debug)]
#[clap(
    name = "bed2gtf",
    version = "1.9.1",
    author = "Alejandro Gonzales-Irribarren <jose.gonzalesdezavala1@unmsm.edu.pe>",
    about = "A fast and memory efficient BED to GTF converter"
)]
pub struct Cli {
    #[clap(
        short = 'b',
        long,
        help = "Path to BED file",
        value_name = "BED",
        required = true
    )]
    pub bed: PathBuf,

    #[clap(
        short = 'o',
        long,
        help = "Path to output file",
        value_name = "OUTPUT",
        required = true
    )]
    pub output: PathBuf,

    #[clap(
        short = 't',
        long,
        help = "Number of threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    #[clap(
        short,
        long = "gz",
        help = "Compress output file",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub gz: bool,

    #[arg(
        short,
        long = "no-gene",
        help = "Flag to disable gene_id feature",
        value_name = "FLAG",
        default_missing_value("true"),
        default_value("false"),
        num_args(0..=1),
        require_equals(true),
        action = ArgAction::Set,
    )]
    pub no_gene: bool,

    #[clap(
        short = 'i',
        long,
        help = "Path to isoforms file",
        value_name = "ISOFORMS",
        required_unless_present = "no_gene",
        default_value = None,
    )]
    pub isoforms: Option<PathBuf>,
}

#[derive(Debug, Error)]
pub enum CliError {
    #[error("Invalid input: {0}")]
    InvalidInput(String),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

impl Cli {
    pub fn check(&self) -> Result<(), CliError> {
        self.validate_args()
    }

    fn validate_args(&self) -> Result<(), CliError> {
        validate(&self.bed)?;

        match self.bed.extension() {
            Some(ext) if ext == "bed" => (),
            _ => {
                return Err(CliError::InvalidInput(format!(
                    "file {:?} is not a BED file",
                    self.bed
                )))
            }
        }

        if !self.no_gene {
            let isoforms = self.isoforms.as_ref().unwrap();
            validate(isoforms)?;
        }

        match self.output.extension() {
            Some(ext) if ext == "gtf" => (),
            _ => {
                return Err(CliError::InvalidInput(format!(
                    "file {:?} is not a GTF file",
                    self.bed
                )))
            }
        }

        Ok(())
    }
}

fn validate(arg: &PathBuf) -> Result<(), CliError> {
    if !arg.exists() {
        return Err(CliError::InvalidInput(format!("{:?} does not exist", arg)));
    }

    if !arg.is_file() {
        return Err(CliError::InvalidInput(format!("{:?} is not a file", arg)));
    }

    match std::fs::metadata(arg) {
        Ok(metadata) if metadata.len() == 0 => {
            Err(CliError::InvalidInput(format!("file {:?} is empty", arg)))
        }
        Ok(_) => Ok(()),
        Err(e) => Err(CliError::IoError(e)),
    }
}
