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
    let mut outputs = reader
        .par_chunks(config.chunks)?
        .map(|(idx, chunk)| render_chunk::<K>(idx, chunk, isoforms))
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

/// Renders a chunk of GenePred records into GXF lines.
///
/// ```rust,ignore
/// # use std::collections::HashMap;
/// # let isoforms: Option<&HashMap<String, String>> = None;
/// # let chunk: Vec<genepred::ReaderResult<genepred::GenePred>> = vec![];
/// let (idx, buffer) = bed2gtf::convert::render_chunk::<genepred::Gtf>(0, chunk, isoforms)?;
/// ```
fn render_chunk<K>(
    idx: usize,
    chunk: Vec<genepred::ReaderResult<GenePred>>,
    isoforms: Option<&HashMap<String, String>>,
) -> Result<(usize, Vec<u8>)>
where
    K: genepred::BedFormat,
{
    let mut buffer = Vec::with_capacity(chunk.len().saturating_mul(256));
    {
        let mut writer = BufWriter::with_capacity(128 * 1024, &mut buffer);
        for record in chunk {
            let record = record?;
            for line in record.to_gxf::<K>(isoforms) {
                writer.write_all(&line)?;
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
