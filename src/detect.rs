// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use crate::cli::OutputFormat;
use crate::error::{Bed2GtfError, Result};
use std::path::{Path, PathBuf};

/// Compression mode used by input and output paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Compression {
    None,
    Gzip,
}

impl Compression {
    /// Returns `true` when the stream is gzip-compressed.
    ///
    /// ```rust,ignore
    /// use bed2gtf::detect::Compression;
    /// assert!(Compression::Gzip.is_gzip());
    /// ```
    pub(crate) fn is_gzip(self) -> bool {
        matches!(self, Compression::Gzip)
    }
}

/// Parsed metadata for a BED input path.
///
/// ```rust,ignore
/// use bed2gtf::detect::{InputPathKind, Compression};
/// use std::path::PathBuf;
/// let kind = InputPathKind { path: PathBuf::from("input.bed"), compression: Compression::None };
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct InputPathKind {
    pub(crate) path: PathBuf,
    pub(crate) compression: Compression,
}

/// Parsed metadata for an output path.
///
/// ```rust,ignore
/// use bed2gtf::detect::{OutputPathKind, Compression};
/// use bed2gtf::OutputFormat;
/// use std::path::PathBuf;
/// let kind = OutputPathKind { path: PathBuf::from("output.gtf"), format: OutputFormat::Gtf, compression: Compression::None };
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct OutputPathKind {
    pub(crate) path: PathBuf,
    pub(crate) format: OutputFormat,
    pub(crate) compression: Compression,
}

/// Detects BED input metadata from a filesystem path.
///
/// ```rust,ignore
/// use bed2gtf::detect::{detect_input_path, Compression};
/// let kind = detect_input_path("input.bed.gz".as_ref())?;
/// assert_eq!(kind.compression, Compression::Gzip);
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
pub(crate) fn detect_input_path(path: &Path) -> Result<InputPathKind> {
    let ext = extension_lowercase(path)
        .ok_or_else(|| Bed2GtfError::UnsupportedInputExtension(path.display().to_string()))?;

    let compression = compression_from_extension(&ext);
    let bed_ext = if compression.is_gzip() {
        nested_extension(path)
            .ok_or_else(|| Bed2GtfError::UnsupportedInputExtension(path.display().to_string()))?
    } else {
        ext
    };

    if bed_ext != "bed" {
        return Err(Bed2GtfError::UnsupportedInputExtension(
            path.display().to_string(),
        ));
    }

    Ok(InputPathKind {
        path: path.to_path_buf(),
        compression,
    })
}

/// Detects output metadata from a filesystem path.
///
/// ```rust,ignore
/// use bed2gtf::{detect::detect_output_path, OutputFormat};
/// let kind = detect_output_path("output.gff3.gz".as_ref())?;
/// assert_eq!(kind.format, OutputFormat::Gff);
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
pub(crate) fn detect_output_path(path: &Path) -> Result<OutputPathKind> {
    let ext = extension_lowercase(path)
        .ok_or_else(|| Bed2GtfError::UnsupportedOutputExtension(path.display().to_string()))?;

    let compression = compression_from_extension(&ext);
    let format_ext = if compression.is_gzip() {
        nested_extension(path)
            .ok_or_else(|| Bed2GtfError::UnsupportedOutputExtension(path.display().to_string()))?
    } else {
        ext
    };

    let format = match format_ext.as_str() {
        "gtf" => OutputFormat::Gtf,
        "gff" | "gff3" => OutputFormat::Gff,
        _ => {
            return Err(Bed2GtfError::UnsupportedOutputExtension(
                path.display().to_string(),
            ))
        }
    };

    Ok(OutputPathKind {
        path: path.to_path_buf(),
        format,
        compression,
    })
}

/// Extracts the lowercase file extension from a path.
///
/// ```rust,ignore
/// use std::path::Path;
/// let ext = bed2gtf::detect::extension_lowercase(Path::new("file.BED"));
/// assert_eq!(ext, Some("bed".to_string()));
/// ```
fn extension_lowercase(path: &Path) -> Option<String> {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.to_ascii_lowercase())
}

/// Extracts the nested extension from a double extension path (e.g., `.bed.gz`).
///
/// ```rust,ignore
/// use std::path::Path;
/// let ext = bed2gtf::detect::nested_extension(Path::new("file.bed.gz"));
/// assert_eq!(ext, Some("bed".to_string()));
/// ```
fn nested_extension(path: &Path) -> Option<String> {
    let stem = path.file_stem()?.to_str()?;
    extension_lowercase(Path::new(stem))
}

/// Determines compression mode from a file extension.
///
/// ```rust,ignore
/// use bed2gtf::detect::{compression_from_extension, Compression};
/// assert_eq!(compression_from_extension("gz"), Compression::Gzip);
/// assert_eq!(compression_from_extension("bed"), Compression::None);
/// ```
fn compression_from_extension(ext: &str) -> Compression {
    match ext {
        "gz" | "gzip" => Compression::Gzip,
        _ => Compression::None,
    }
}

#[cfg(test)]
mod tests {
    use super::{detect_input_path, detect_output_path, Compression};
    use crate::OutputFormat;
    use std::path::Path;

    #[test]
    fn detects_bed_gzip_input() {
        let kind = detect_input_path(Path::new("input.bed.gz")).unwrap();
        assert_eq!(kind.compression, Compression::Gzip);
    }

    #[test]
    fn detects_gff3_output() {
        let kind = detect_output_path(Path::new("out.gff3.gz")).unwrap();
        assert_eq!(kind.format, OutputFormat::Gff);
        assert_eq!(kind.compression, Compression::Gzip);
    }
}
