// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use crate::cli::IsoformFormat;
use crate::error::{Bed2GtfError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Loads an optional transcript-to-gene map from disk.
///
/// ```rust,ignore
/// use bed2gtf::isoforms::load_isoforms;
/// let mapping = load_isoforms(None, bed2gtf::cli::IsoformFormat::Gt)?;
/// assert!(mapping.is_none());
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
pub(crate) fn load_isoforms(
    path: Option<&Path>,
    format: IsoformFormat,
) -> Result<Option<HashMap<String, String>>> {
    path.map(|p| parse_isoforms(p, format)).transpose()
}

/// Parses a two-column tab-separated transcript-to-gene file.
///
/// The default format is `g,t` (gene first, as written by TOGA2); pass
/// `IsoformFormat::Tg` for transcript-first files.
///
/// ```rust,ignore
/// use bed2gtf::isoforms::parse_isoforms;
/// let mapping = parse_isoforms("isoforms.tsv".as_ref(), bed2gtf::cli::IsoformFormat::Gt)?;
/// assert!(mapping.contains_key("ENST0001"));
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
pub(crate) fn parse_isoforms(
    path: &Path,
    format: IsoformFormat,
) -> Result<HashMap<String, String>> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(128 * 1024, file);
    let mut mapping = HashMap::new();

    for (index, line) in reader.lines().enumerate() {
        let line_no = index + 1;
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            return Err(Bed2GtfError::InvalidIsoformsLine {
                path: path.to_path_buf(),
                line: line_no,
                message: "expected exactly two tab-separated columns".into(),
            });
        }

        let (transcript, gene) = match format {
            // Gene first: `gene<TAB>transcript` (v1.x / TOGA2).
            IsoformFormat::Gt => (fields[1].trim(), fields[0].trim()),
            // Transcript first: `transcript<TAB>gene` (v2.1.0).
            IsoformFormat::Tg => (fields[0].trim(), fields[1].trim()),
        };
        if transcript.is_empty() || gene.is_empty() {
            return Err(Bed2GtfError::InvalidIsoformsLine {
                path: path.to_path_buf(),
                line: line_no,
                message: "transcript and gene identifiers must be non-empty".into(),
            });
        }

        if mapping
            .insert(transcript.to_owned(), gene.to_owned())
            .is_some()
        {
            return Err(Bed2GtfError::DuplicateIsoformMapping {
                path: path.to_path_buf(),
                line: line_no,
                transcript: transcript.to_owned(),
            });
        }
    }

    Ok(mapping)
}

#[cfg(test)]
mod tests {
    use super::parse_isoforms;
    use crate::cli::IsoformFormat;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(name: &str) -> std::path::PathBuf {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("bed2gtf-{name}-{stamp}.tsv"))
    }

    #[test]
    fn parses_gene_first_map() {
        let path = temp_path("ok");
        fs::write(&path, "gene1\ttx1\ngene2\ttx2\n").unwrap();

        let mapping = parse_isoforms(&path, IsoformFormat::Gt).unwrap();
        assert_eq!(mapping.get("tx1").unwrap(), "gene1");

        let _ = fs::remove_file(path);
    }

    #[test]
    fn parses_transcript_first_map() {
        let path = temp_path("tg");
        fs::write(&path, "tx1\tgene1\ntx2\tgene2\n").unwrap();

        let mapping = parse_isoforms(&path, IsoformFormat::Tg).unwrap();
        assert_eq!(mapping.get("tx1").unwrap(), "gene1");

        let _ = fs::remove_file(path);
    }

    #[test]
    fn rejects_duplicate_transcripts() {
        let path = temp_path("dup");
        fs::write(&path, "gene1\ttx1\ngene2\ttx1\n").unwrap();

        assert!(parse_isoforms(&path, IsoformFormat::Gt).is_err());

        let _ = fs::remove_file(path);
    }
}
