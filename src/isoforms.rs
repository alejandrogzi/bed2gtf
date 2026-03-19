// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use crate::error::{Bed2GtfError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Loads an optional transcript-to-gene map from disk.
///
/// ```rust,ignore
/// use bed2gtf::isoforms::load_isoforms;
/// let mapping = load_isoforms(None)?;
/// assert!(mapping.is_none());
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
pub(crate) fn load_isoforms(path: Option<&Path>) -> Result<Option<HashMap<String, String>>> {
    path.map(parse_isoforms).transpose()
}

/// Parses a two-column tab-separated transcript-to-gene file.
///
/// ```rust,ignore
/// use bed2gtf::isoforms::parse_isoforms;
/// let mapping = parse_isoforms("isoforms.tsv".as_ref())?;
/// assert!(mapping.contains_key("ENST0001"));
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
pub(crate) fn parse_isoforms(path: &Path) -> Result<HashMap<String, String>> {
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

        let transcript = fields[0].trim();
        let gene = fields[1].trim();
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
    fn parses_isoform_map() {
        let path = temp_path("ok");
        fs::write(&path, "tx1\tgene1\ntx2\tgene2\n").unwrap();

        let mapping = parse_isoforms(&path).unwrap();
        assert_eq!(mapping.get("tx1").unwrap(), "gene1");

        let _ = fs::remove_file(path);
    }

    #[test]
    fn rejects_duplicate_transcripts() {
        let path = temp_path("dup");
        fs::write(&path, "tx1\tgene1\ntx1\tgene2\n").unwrap();

        assert!(parse_isoforms(&path).is_err());

        let _ = fs::remove_file(path);
    }
}
