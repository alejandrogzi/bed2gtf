// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::{ArgAction, Parser, ValueEnum};
use log::Level;
use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Parser, Debug, Clone)]
#[command(
    name = "bed2gtf",
    version = env!("CARGO_PKG_VERSION"),
    author = "Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>",
    about = "A fast and memory efficient BED to GTF/GFF converter"
)]
pub struct Args {
    /// Input BED path. Reads from stdin when omitted.
    #[arg(short = 'i', long = "input", value_name = "BED")]
    pub input: Option<PathBuf>,

    /// Output path. Writes to stdout when omitted.
    #[arg(short = 'o', long = "output", value_name = "OUTPUT")]
    pub output: Option<PathBuf>,

    /// Output format when writing to stdout.
    #[arg(long = "to", value_name = "FORMAT")]
    pub to: Option<OutputFormat>,

    /// Optional transcript-to-gene mapping file.
    #[arg(short = 'I', long = "isoforms", value_name = "TSV")]
    pub isoforms: Option<PathBuf>,

    /// BED layout to parse.
    #[arg(
        short = 't',
        long = "type",
        value_name = "BED_TYPE",
        default_value_t = BedType::Bed12
    )]
    pub bed_type: BedType,

    /// Number of worker threads.
    #[arg(
        short = 'T',
        long = "threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get()
    )]
    pub threads: usize,

    /// Chunk size for parallel processing.
    #[arg(
        short = 'c',
        long = "chunks",
        value_name = "CHUNKS",
        default_value_t = 15_000
    )]
    pub chunks: usize,

    /// Gzip stdout or require a gzip output path.
    #[arg(short = 'g', long = "gz", action = ArgAction::SetTrue)]
    pub gz: bool,

    /// Logging verbosity.
    #[arg(
        short = 'L',
        long = "level",
        value_name = "LEVEL",
        default_value_t = Level::Info
    )]
    pub level: Level,
}

/// Output layouts supported by the converter.
///
/// ```rust,ignore
/// use bed2gtf::OutputFormat;
/// use std::str::FromStr;
///
/// assert_eq!(OutputFormat::from_str("gff").unwrap(), OutputFormat::Gff);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum OutputFormat {
    Gtf,
    Gff,
}

impl fmt::Display for OutputFormat {
    /// Returns the CLI name of the output format.
    ///
    /// ```rust,ignore
    /// use bed2gtf::OutputFormat;
    /// assert_eq!(OutputFormat::Gtf.to_string(), "gtf");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OutputFormat::Gtf => f.write_str("gtf"),
            OutputFormat::Gff => f.write_str("gff"),
        }
    }
}

/// BED layouts supported by the reader.
///
/// ```rust,ignore
/// use bed2gtf::BedType;
/// assert_eq!(BedType::default(), BedType::Bed12);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BedType {
    Bed3,
    Bed4,
    Bed5,
    Bed6,
    Bed8,
    Bed9,
    Bed12,
}

impl Default for BedType {
    /// Returns the default BED layout.
    ///
    /// ```rust,ignore
    /// use bed2gtf::BedType;
    /// assert_eq!(BedType::default(), BedType::Bed12);
    /// ```
    fn default() -> Self {
        BedType::Bed12
    }
}

impl FromStr for BedType {
    type Err = String;

    /// Parses a BED layout from its numeric representation.
    ///
    /// ```rust,ignore
    /// use bed2gtf::BedType;
    /// use std::str::FromStr;
    ///
    /// assert_eq!(BedType::from_str("8").unwrap(), BedType::Bed8);
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "3" => Ok(BedType::Bed3),
            "4" => Ok(BedType::Bed4),
            "5" => Ok(BedType::Bed5),
            "6" => Ok(BedType::Bed6),
            "8" => Ok(BedType::Bed8),
            "9" => Ok(BedType::Bed9),
            "12" => Ok(BedType::Bed12),
            _ => Err(format!("invalid BED type: {s}")),
        }
    }
}

impl fmt::Display for BedType {
    /// Formats the BED layout as its numeric form.
    ///
    /// ```rust,ignore
    /// use bed2gtf::BedType;
    /// assert_eq!(BedType::Bed6.to_string(), "6");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BedType::Bed3 => f.write_str("3"),
            BedType::Bed4 => f.write_str("4"),
            BedType::Bed5 => f.write_str("5"),
            BedType::Bed6 => f.write_str("6"),
            BedType::Bed8 => f.write_str("8"),
            BedType::Bed9 => f.write_str("9"),
            BedType::Bed12 => f.write_str("12"),
        }
    }
}
