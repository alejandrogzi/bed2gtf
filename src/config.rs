// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use crate::cli::{Args, BedType, OutputFormat};
use crate::detect::{
    detect_input_path, detect_output_path, Compression, InputPathKind, OutputPathKind,
};
use crate::error::{Bed2GtfError, Result};
use std::path::PathBuf;

/// Normalized runtime configuration.
#[derive(Clone, Debug)]
pub(crate) struct Config {
    pub(crate) input: InputSource,
    pub(crate) output: OutputTarget,
    pub(crate) output_format: OutputFormat,
    pub(crate) bed_type: BedType,
    pub(crate) threads: usize,
    pub(crate) chunks: usize,
    pub(crate) isoforms: Option<PathBuf>,
    pub(crate) force_gzip: bool,
}

/// Input source used by a conversion run.
#[derive(Clone, Debug)]
pub(crate) enum InputSource {
    Path(InputPath),
    Stdin,
}

/// Parsed BED input path metadata.
#[derive(Clone, Debug)]
pub(crate) struct InputPath {
    pub(crate) path: PathBuf,
    pub(crate) compression: Compression,
}

/// Output destination used by a conversion run.
#[derive(Clone, Debug)]
pub(crate) enum OutputTarget {
    Path(OutputPath),
    Stdout,
}

/// Parsed output path metadata.
#[derive(Clone, Debug)]
pub(crate) struct OutputPath {
    pub(crate) path: PathBuf,
    pub(crate) compression: Compression,
}

impl Config {
    /// Builds a normalized config from CLI arguments.
    ///
    /// ```rust,ignore
    /// use bed2gtf::{Args, BedType};
    /// let args = Args::parse_from(["bed2gtf", "--to", "gtf"]);
    /// let config = bed2gtf::config::Config::from_args(args)?;
    /// assert_eq!(config.bed_type, BedType::Bed12);
    /// # Ok::<(), bed2gtf::Bed2GtfError>(())
    /// ```
    pub(crate) fn from_args(args: Args) -> Result<Self> {
        if args.threads == 0 {
            return Err(Bed2GtfError::InvalidThreads(args.threads));
        }
        if args.chunks == 0 {
            return Err(Bed2GtfError::InvalidChunkSize(args.chunks));
        }

        let input = match args.input {
            Some(path) => {
                let kind = detect_input_path(&path)?;
                InputSource::Path(InputPath::from(kind))
            }
            None => InputSource::Stdin,
        };

        let (output, output_format) = match args.output {
            Some(path) => {
                if args.to.is_some() {
                    return Err(Bed2GtfError::ConflictingOutputFormat);
                }

                let kind = detect_output_path(&path)?;
                validate_output_path(&input, &kind, args.gz)?;
                let format = kind.format;
                (OutputTarget::Path(OutputPath::from(kind)), format)
            }
            None => {
                let format = args.to.ok_or(Bed2GtfError::MissingOutputFormat)?;
                (OutputTarget::Stdout, format)
            }
        };

        Ok(Self {
            input,
            output,
            output_format,
            bed_type: args.bed_type,
            threads: args.threads,
            chunks: args.chunks,
            isoforms: args.isoforms,
            force_gzip: args.gz,
        })
    }
}

impl From<InputPathKind> for InputPath {
    /// Converts `InputPathKind` to `InputPath` by extracting path and compression.
    ///
    /// ```rust,ignore
    /// # use bed2gtf::detect::InputPathKind;
    /// # use std::path::PathBuf;
    /// # use bed2gtf::detect::Compression;
    /// let kind = InputPathKind { path: PathBuf::from("input.bed"), compression: Compression::None };
    /// let input: bed2gtf::config::InputPath = kind.into();
    /// ```
    fn from(kind: InputPathKind) -> Self {
        Self {
            path: kind.path,
            compression: kind.compression,
        }
    }
}

impl From<OutputPathKind> for OutputPath {
    /// Converts `OutputPathKind` to `OutputPath` by extracting path and compression.
    ///
    /// ```rust,ignore
    /// # use bed2gtf::detect::OutputPathKind;
    /// # use std::path::PathBuf;
    /// # use bed2gtf::detect::Compression;
    /// # use bed2gtf::OutputFormat;
    /// let kind = OutputPathKind { path: PathBuf::from("output.gtf"), format: OutputFormat::Gtf, compression: Compression::None };
    /// let output: bed2gtf::config::OutputPath = kind.into();
    /// ```
    fn from(kind: OutputPathKind) -> Self {
        Self {
            path: kind.path,
            compression: kind.compression,
        }
    }
}

/// Validates output path settings against input source.
///
/// Ensures gzip flag consistency and that compressed input requires gzip output.
///
/// ```rust,ignore
/// # use bed2gtf::config::{InputSource, OutputTarget, OutputPath};
/// # use bed2gtf::detect::{OutputPathKind, Compression};
/// # use bed2gtf::OutputFormat;
/// # use std::path::PathBuf;
/// let input = InputSource::Stdin;
/// let output = OutputPathKind { path: PathBuf::from("out.gtf"), format: OutputFormat::Gtf, compression: Compression::None };
/// bed2gtf::config::validate_output_path(&input, &output, false)?;
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
fn validate_output_path(input: &InputSource, output: &OutputPathKind, gz_flag: bool) -> Result<()> {
    if gz_flag && !output.compression.is_gzip() {
        return Err(Bed2GtfError::OutputPathMustBeGzip(output.path.clone()));
    }

    if let InputSource::Path(input) = input {
        if input.compression.is_gzip() && !output.compression.is_gzip() {
            return Err(Bed2GtfError::CompressedInputRequiresGzipOutput {
                input: input.path.clone(),
                output: output.path.clone(),
            });
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{Config, InputSource, OutputTarget};
    use crate::{Args, OutputFormat};
    use clap::Parser;

    #[test]
    fn requires_stdout_format_without_output_path() {
        let args = Args::parse_from(["bed2gtf"]);
        let err = Config::from_args(args).unwrap_err();
        assert_eq!(
            err.to_string(),
            "missing output format: use --output or --to when writing to stdout"
        );
    }

    #[test]
    fn resolves_stdout_target() {
        let args = Args::parse_from(["bed2gtf", "--to", "gff"]);
        let config = Config::from_args(args).unwrap();
        assert!(matches!(config.input, InputSource::Stdin));
        assert!(matches!(config.output, OutputTarget::Stdout));
        assert_eq!(config.output_format, OutputFormat::Gff);
    }
}
