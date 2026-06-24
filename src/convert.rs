// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use crate::cli::{BedType, OutputFormat};
use crate::config::{Config, InputSource, OutputTarget};
use crate::detect::Compression;
use crate::error::Result;
use crate::isoforms::load_isoforms;
use crate::memory::max_mem_usage_mb;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression as GzipLevel;
use genepred::{Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, GenePred, Gff, Gtf, Reader};
use rayon::prelude::*;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Cursor, Read, Write};
use std::time::{Duration, Instant};

/// Summary statistics for a conversion run.
#[derive(Debug, Clone, Copy)]
pub(crate) struct RunStats {
    pub(crate) elapsed: Duration,
    pub(crate) mem_delta_mb: f64,
}

/// Runs one BED-to-GTF/GFF conversion.
///
/// ```rust,ignore
/// # use bed2gtf::config::Config;
/// # let config: Config = todo!();
/// let _stats = bed2gtf::convert::run(&config)?;
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
pub(crate) fn run(config: &Config) -> Result<RunStats> {
    let start = Instant::now();
    let start_mem = max_mem_usage_mb();
    let isoforms = load_isoforms(config.isoforms.as_deref())?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build()?;

    let (input_compression, chunks) = match config.output_format {
        OutputFormat::Gtf => process_by_bed_type::<Gtf>(config, isoforms.as_ref(), &pool)?,
        OutputFormat::Gff => process_by_bed_type::<Gff>(config, isoforms.as_ref(), &pool)?,
    };

    validate_runtime_output(config, input_compression)?;
    let sink = open_output_sink(config, input_compression)?;
    let mut writer = BufWriter::with_capacity(256 * 1024, sink);
    write_output(&mut writer, chunks)?;

    Ok(RunStats {
        elapsed: start.elapsed(),
        mem_delta_mb: (max_mem_usage_mb() - start_mem).max(0.0),
    })
}

/// Dispatches conversion to the correct BED format handler.
///
/// ```rust,ignore
/// # use bed2gtf::config::Config;
/// # let config: Config = todo!();
/// # let isoforms: Option<&std::collections::HashMap<String, String>> = None;
/// # let pool: rayon::ThreadPool = todo!();
/// let _ = bed2gtf::convert::process_by_bed_type::<bed2gtf::cli::OutputFormat>(
///     &config,
///     isoforms,
///     &pool,
/// );
/// ```
#[allow(clippy::type_complexity)]
fn process_by_bed_type<K>(
    config: &Config,
    isoforms: Option<&HashMap<String, String>>,
    pool: &rayon::ThreadPool,
) -> Result<(Compression, Vec<(usize, Vec<u8>)>)>
where
    K: genepred::BedFormat,
{
    match config.bed_type {
        BedType::Bed3 => pool.install(|| process::<Bed3, K>(config, isoforms)),
        BedType::Bed4 => pool.install(|| process::<Bed4, K>(config, isoforms)),
        BedType::Bed5 => pool.install(|| process::<Bed5, K>(config, isoforms)),
        BedType::Bed6 => pool.install(|| process::<Bed6, K>(config, isoforms)),
        BedType::Bed8 => pool.install(|| process::<Bed8, K>(config, isoforms)),
        BedType::Bed9 => pool.install(|| process::<Bed9, K>(config, isoforms)),
        BedType::Bed12 => pool.install(|| process::<Bed12, K>(config, isoforms)),
    }
}

/// Converts a BED chunk into GXF records using parallel processing.
///
/// ```rust,ignore
/// # use bed2gtf::config::Config;
/// # let config: Config = todo!();
/// # let isoforms: Option<&std::collections::HashMap<String, String>> = None;
/// let (compression, chunks) = bed2gtf::convert::process::<genepred::Bed12, genepred::Gtf>(
///     &config,
///     isoforms,
/// )?;
/// ```
#[allow(clippy::type_complexity)]
fn process<B, K>(
    config: &Config,
    isoforms: Option<&HashMap<String, String>>,
) -> Result<(Compression, Vec<(usize, Vec<u8>)>)>
where
    B: genepred::BedFormat + Into<GenePred> + Send,
    K: genepred::BedFormat,
{
    let (reader, input_compression) = open_reader::<B>(config)?;

    // Materialize every record in input order: parse chunks in parallel, then
    // order them by chunk index before flattening.
    let mut chunks = reader
        .par_chunks(config.chunks)?
        .collect::<Vec<(usize, Vec<genepred::ReaderResult<GenePred>>)>>();
    chunks.sort_by_key(|(idx, _)| *idx);

    let mut records = Vec::new();
    for (_, chunk) in chunks {
        for record in chunk {
            records.push(record?);
        }
    }

    // First pass: aggregate each gene's span (minimum start, maximum end) across
    // its isoforms and flag, for every record, whether it is the first
    // occurrence of its gene (its gene "leader"). genepred emits a `gene` line
    // per record, so without this only the leader keeps a single, span-correct
    // gene line; the others have their duplicate gene line dropped during render.
    let mut spans: HashMap<Vec<u8>, (u64, u64)> = HashMap::new();
    let mut leaders = Vec::with_capacity(records.len());
    for record in &records {
        match resolve_gene_key(record, isoforms) {
            Some(key) => match spans.entry(key) {
                Entry::Occupied(mut slot) => {
                    let (start, end) = slot.get_mut();
                    *start = (*start).min(record.start);
                    *end = (*end).max(record.end);
                    leaders.push(false);
                }
                Entry::Vacant(slot) => {
                    slot.insert((record.start, record.end));
                    leaders.push(true);
                }
            },
            None => leaders.push(true),
        }
    }

    // Second pass: render windows of records in parallel, emitting exactly one
    // span-corrected gene line per gene at its leader's position.
    let window = config.chunks.max(1);
    let mut outputs = (0..records.len())
        .step_by(window)
        .enumerate()
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|(idx, base)| {
            let end = (base + window).min(records.len());
            render_window::<K>(idx, base, &records[base..end], &leaders, &spans, isoforms)
        })
        .collect::<Vec<_>>();

    let mut merged = Vec::with_capacity(outputs.len());
    for chunk in outputs.drain(..) {
        merged.push(chunk?);
    }

    merged.sort_by_key(|(idx, _)| *idx);
    Ok((input_compression, merged))
}

/// Opens a BED reader from a file path or stdin.
///
/// ```rust,ignore
/// # use bed2gtf::config::Config;
/// # let config: Config = todo!();
/// let (reader, compression) = bed2gtf::convert::open_reader::<genepred::Bed12>(&config)?;
/// ```
fn open_reader<B>(config: &Config) -> Result<(Reader<B>, Compression)>
where
    B: genepred::BedFormat + Into<GenePred> + Send,
{
    match &config.input {
        InputSource::Path(input) => {
            let reader = if input.compression.is_gzip() {
                Reader::<B>::from_path(&input.path)?
            } else {
                Reader::<B>::from_mmap(&input.path)?
            };
            Ok((reader, input.compression))
        }
        InputSource::Stdin => {
            let (reader, compression) = open_stdin_reader()?;
            let reader = Reader::<B>::from_reader(reader)?;
            Ok((reader, compression))
        }
    }
}

/// Opens a stdin reader with automatic gzip detection.
///
/// ```rust,ignore
/// // Reads first 2 bytes to detect gzip magic header
/// let (reader, compression) = bed2gtf::convert::open_stdin_reader()?;
/// ```
fn open_stdin_reader() -> Result<(Box<dyn Read + Send>, Compression)> {
    let mut stdin = io::stdin();
    let mut prefix = [0u8; 2];
    let bytes = stdin.read(&mut prefix)?;
    let prefix_reader = Cursor::new(prefix[..bytes].to_vec());
    let chained = prefix_reader.chain(stdin);

    if bytes == 2 && prefix == [0x1f, 0x8b] {
        Ok((Box::new(MultiGzDecoder::new(chained)), Compression::Gzip))
    } else {
        Ok((Box::new(chained), Compression::None))
    }
}

/// Resolves the gene grouping key for a record, mirroring how genepred resolves
/// the gene identifier it stamps on the emitted GXF lines.
///
/// Returns the mapped gene when an isoform mapping is supplied and matches the
/// record name, the record name otherwise, or `None` when the record has no
/// name (in which case it cannot be grouped and is treated as its own gene).
///
/// ```rust,ignore
/// # use std::collections::HashMap;
/// # let record: genepred::GenePred = todo!();
/// let key = bed2gtf::convert::resolve_gene_key(&record, None);
/// ```
fn resolve_gene_key(
    record: &GenePred,
    isoforms: Option<&HashMap<String, String>>,
) -> Option<Vec<u8>> {
    let name = record.name.as_ref()?;
    if let Some(mapping) = isoforms {
        if let Ok(name_text) = std::str::from_utf8(name) {
            if let Some(gene) = mapping.get(name_text) {
                return Some(gene.as_bytes().to_vec());
            }
        }
    }
    Some(name.clone())
}

/// Returns true when a GXF line's feature column is `gene`.
fn gene_feature_is_gene(line: &[u8]) -> bool {
    line.split(|&byte| byte == b'\t').nth(2) == Some(b"gene".as_ref())
}

/// Rewrites the start and end coordinate columns of a GXF `gene` line.
///
/// The feature spans an entire gene, so its coordinates must reflect the union
/// of all isoforms rather than the single record genepred derived it from.
///
/// ```rust,ignore
/// let line = b"chr1\tTOGA2\tgene\t11\t20\t.\t+\t.\tgene_id \"g\";";
/// let fixed = bed2gtf::convert::rewrite_gene_span(line, 5, 40);
/// ```
fn rewrite_gene_span(line: &[u8], start_col: u64, end_col: u64) -> Vec<u8> {
    let start = start_col.to_string();
    let end = end_col.to_string();
    let mut out = Vec::with_capacity(line.len() + start.len() + end.len());
    for (index, field) in line.split(|&byte| byte == b'\t').enumerate() {
        if index > 0 {
            out.push(b'\t');
        }
        match index {
            3 => out.extend_from_slice(start.as_bytes()),
            4 => out.extend_from_slice(end.as_bytes()),
            _ => out.extend_from_slice(field),
        }
    }
    out
}

/// Renders a window of GenePred records into GXF lines.
///
/// genepred emits a `gene` line for every record; this drops that per-record
/// gene line and emits a single span-corrected gene line only when the record
/// is its gene's leader (`leaders[base + offset]`).
///
/// ```rust,ignore
/// # use std::collections::HashMap;
/// # let isoforms: Option<&HashMap<String, String>> = None;
/// # let records: Vec<genepred::GenePred> = vec![];
/// # let leaders = vec![true];
/// # let spans = HashMap::new();
/// let (idx, buffer) =
///     bed2gtf::convert::render_window::<genepred::Gtf>(0, 0, &records, &leaders, &spans, isoforms)?;
/// ```
fn render_window<K>(
    idx: usize,
    base: usize,
    records: &[GenePred],
    leaders: &[bool],
    spans: &HashMap<Vec<u8>, (u64, u64)>,
    isoforms: Option<&HashMap<String, String>>,
) -> Result<(usize, Vec<u8>)>
where
    K: genepred::BedFormat,
{
    let mut buffer = Vec::with_capacity(records.len().saturating_mul(256));
    {
        let mut writer = BufWriter::with_capacity(128 * 1024, &mut buffer);
        for (offset, record) in records.iter().enumerate() {
            let lines = record.to_gxf::<K>(isoforms);
            let mut lines = lines.iter();

            // genepred always emits the per-record `gene` line first; drop it so
            // a gene appears exactly once, contributed by its leader below.
            let gene_line = lines.next();
            debug_assert!(
                gene_line
                    .map(|line| gene_feature_is_gene(line))
                    .unwrap_or(false),
                "expected genepred to emit the gene feature as the first GXF line",
            );

            if leaders[base + offset] {
                if let Some(gene_line) = gene_line {
                    let (start, end) = resolve_gene_key(record, isoforms)
                        .and_then(|key| spans.get(&key).copied())
                        .unwrap_or((record.start, record.end));
                    let corrected = rewrite_gene_span(gene_line, start.saturating_add(1), end);
                    writer.write_all(&corrected)?;
                    writer.write_all(b"\n")?;
                }
            }

            for line in lines {
                writer.write_all(line)?;
                writer.write_all(b"\n")?;
            }
        }
        writer.flush()?;
    }

    Ok((idx, buffer))
}

/// Opens the output sink (file or stdout) with optional gzip compression.
///
/// ```rust,ignore
/// # use bed2gtf::config::Config;
/// # let config: Config = todo!();
/// # use bed2gtf::detect::Compression;
/// let sink = bed2gtf::convert::open_output_sink(&config, Compression::None)?;
/// ```
fn open_output_sink(config: &Config, input_compression: Compression) -> Result<Box<dyn Write>> {
    match &config.output {
        OutputTarget::Path(output) => {
            let file = File::create(&output.path)?;
            if output.compression.is_gzip() {
                Ok(Box::new(GzEncoder::new(file, GzipLevel::fast())))
            } else {
                Ok(Box::new(file))
            }
        }
        OutputTarget::Stdout => {
            let gzip = config.force_gzip || input_compression.is_gzip();
            let stdout = io::stdout();
            if gzip {
                Ok(Box::new(GzEncoder::new(stdout, GzipLevel::fast())))
            } else {
                Ok(Box::new(stdout))
            }
        }
    }
}

/// Writes sorted chunks to the output writer in order.
///
/// ```rust,ignore
/// # use std::io::stdout;
/// # let chunks: Vec<(usize, Vec<u8>)> = vec![];
/// # let mut writer = stdout();
/// bed2gtf::convert::write_output(&mut writer, chunks)?;
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
fn write_output<W: Write>(writer: &mut W, chunks: Vec<(usize, Vec<u8>)>) -> Result<()> {
    for (_, buffer) in chunks {
        writer.write_all(&buffer)?;
    }
    writer.flush()?;
    Ok(())
}

/// Validates that runtime output settings are compatible with input compression.
///
/// Ensures gzip stdin is not piped to non-gzip file output.
///
/// ```rust,ignore
/// # use bed2gtf::config::{Config, InputSource, OutputTarget};
/// # let config = Config { input: InputSource::Stdin, output: OutputTarget::Stdout, output_format: bed2gtf::OutputFormat::Gtf, bed_type: bed2gtf::BedType::Bed12, threads: 1, chunks: 1, isoforms: None, force_gzip: false };
/// # use bed2gtf::detect::Compression;
/// bed2gtf::convert::validate_runtime_output(&config, Compression::None)?;
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
fn validate_runtime_output(config: &Config, input_compression: Compression) -> Result<()> {
    if matches!(config.input, InputSource::Stdin) && input_compression.is_gzip() {
        if let OutputTarget::Path(output) = &config.output {
            if !output.compression.is_gzip() {
                return Err(
                    crate::error::Bed2GtfError::CompressedStdinRequiresGzipOutput(
                        output.path.clone(),
                    ),
                );
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::run;
    use crate::config::{Config, OutputTarget};
    use crate::Args;
    use clap::Parser;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    /// Three isoforms of one gene at overlapping but distinct spans.
    const BED: &str = "chr1\t1000\t2000\ttxA\t0\t+\t1000\t2000\t0,0,0\t1\t1000,\t0,\n\
chr1\t1100\t2200\ttxB\t0\t+\t1100\t2200\t0,0,0\t1\t1100,\t0,\n\
chr1\t900\t1800\ttxC\t0\t+\t900\t1800\t0,0,0\t1\t900,\t0,\n";

    const ISOFORMS: &str = "txA\tGENE1\ntxB\tGENE1\ntxC\tGENE1\n";

    fn temp_path(tag: &str, ext: &str) -> PathBuf {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("bed2gtf-convert-{tag}-{stamp}.{ext}"))
    }

    /// Runs a full conversion through the CLI surface and returns the output text.
    fn convert(args: &[&str]) -> String {
        let config = Config::from_args(Args::parse_from(args)).unwrap();
        let OutputTarget::Path(output) = &config.output else {
            panic!("test expects a file output target");
        };
        let path = output.path.clone();
        run(&config).unwrap();
        let text = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        text
    }

    fn gene_lines(text: &str) -> Vec<&str> {
        text.lines()
            .filter(|line| line.split('\t').nth(2) == Some("gene"))
            .collect()
    }

    fn write_inputs(tag: &str) -> (PathBuf, PathBuf) {
        let bed = temp_path(tag, "bed");
        let iso = temp_path(tag, "tsv");
        fs::write(&bed, BED).unwrap();
        fs::write(&iso, ISOFORMS).unwrap();
        (bed, iso)
    }

    #[test]
    fn collapses_isoforms_into_single_gene_line_gtf() {
        let (bed, iso) = write_inputs("gtf");
        let out = temp_path("gtf", "gtf");
        let text = convert(&[
            "bed2gtf",
            "-i",
            bed.to_str().unwrap(),
            "-I",
            iso.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ]);
        let _ = fs::remove_file(&bed);
        let _ = fs::remove_file(&iso);

        let genes = gene_lines(&text);
        assert_eq!(genes.len(), 1, "expected one gene line, got: {genes:?}");
        let fields: Vec<&str> = genes[0].split('\t').collect();
        // Union span: 0-based min start 900 -> column 901; 1-based max end 2200.
        assert_eq!(fields[3], "901");
        assert_eq!(fields[4], "2200");
        // The single gene line precedes all transcript rows for that gene.
        let first_feature = text.lines().next().unwrap().split('\t').nth(2).unwrap();
        assert_eq!(first_feature, "gene");
    }

    #[test]
    fn collapses_isoforms_into_single_gene_line_gff() {
        let (bed, iso) = write_inputs("gff");
        let out = temp_path("gff", "gff");
        let text = convert(&[
            "bed2gtf",
            "-i",
            bed.to_str().unwrap(),
            "-I",
            iso.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ]);
        let _ = fs::remove_file(&bed);
        let _ = fs::remove_file(&iso);

        let genes = gene_lines(&text);
        assert_eq!(genes.len(), 1, "expected one gene line, got: {genes:?}");
        let fields: Vec<&str> = genes[0].split('\t').collect();
        assert_eq!(fields[3], "901");
        assert_eq!(fields[4], "2200");
    }

    #[test]
    fn dedup_holds_when_isoforms_span_separate_windows() {
        let (bed, iso) = write_inputs("cross-window");
        let out = temp_path("cross-window", "gtf");
        // chunk size 1 places each isoform in its own parallel window, so a single
        // gene line proves the leader is decided globally, not per window.
        let text = convert(&[
            "bed2gtf",
            "-i",
            bed.to_str().unwrap(),
            "-I",
            iso.to_str().unwrap(),
            "-c",
            "1",
            "-o",
            out.to_str().unwrap(),
        ]);
        let _ = fs::remove_file(&bed);
        let _ = fs::remove_file(&iso);

        let genes = gene_lines(&text);
        assert_eq!(genes.len(), 1, "expected one gene line, got: {genes:?}");
        let fields: Vec<&str> = genes[0].split('\t').collect();
        assert_eq!(fields[3], "901");
        assert_eq!(fields[4], "2200");
    }

    #[test]
    fn keeps_one_gene_line_per_record_without_isoforms() {
        let bed = temp_path("noiso", "bed");
        fs::write(&bed, BED).unwrap();
        let out = temp_path("noiso", "gtf");
        let text = convert(&[
            "bed2gtf",
            "-i",
            bed.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ]);
        let _ = fs::remove_file(&bed);

        // Each transcript is its own gene, so the per-record gene line is correct.
        assert_eq!(gene_lines(&text).len(), 3);
    }

    /// In GFF, an unmapped transcript must not reuse its own id as the gene id:
    /// that collided the `gene`/`mRNA` `ID`s and made the mRNA its own `Parent`
    /// (a self-loop for GFF3 parsers). genepred now synthesizes `gene-<tx>`.
    #[test]
    fn gff_without_isoforms_uses_distinct_gene_and_transcript_ids() {
        let bed = temp_path("noiso-gff", "bed");
        fs::write(&bed, BED).unwrap();
        let out = temp_path("noiso-gff", "gff");
        let text = convert(&[
            "bed2gtf",
            "-i",
            bed.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ]);
        let _ = fs::remove_file(&bed);

        let attr = |line: &str, key: &str| -> Option<String> {
            line.split('\t')
                .nth(8)?
                .split(';')
                .find_map(|a| a.strip_prefix(key).map(str::to_string))
        };

        let mut mrnas = 0;
        for line in text.lines() {
            if line.split('\t').nth(2) == Some("mRNA") {
                mrnas += 1;
                let id = attr(line, "ID=").expect("mRNA has ID");
                let parent = attr(line, "Parent=").expect("mRNA has Parent");
                assert_ne!(parent, id, "mRNA must not be its own parent: {line}");
                assert_eq!(parent, format!("gene-{id}"), "line: {line}");
            }
        }
        assert_eq!(mrnas, 3);

        assert!(gene_lines(&text)
            .iter()
            .filter_map(|line| attr(line, "ID="))
            .all(|gene_id| gene_id.starts_with("gene-")));
    }
}
