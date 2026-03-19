// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use std::path::PathBuf;

use thiserror::Error;

/// Result type used across the crate.
pub type Result<T> = std::result::Result<T, Bed2GtfError>;

/// Error type for BED-to-GTF/GFF conversion.
#[derive(Debug, Error)]
pub enum Bed2GtfError {
    #[error("missing output format: use --output or --to when writing to stdout")]
    MissingOutputFormat,
    #[error("invalid arguments: --to can only be used when --output is absent")]
    ConflictingOutputFormat,
    #[error("invalid thread count: {0}")]
    InvalidThreads(usize),
    #[error("invalid chunk size: {0}")]
    InvalidChunkSize(usize),
    #[error("unsupported input extension: {0}")]
    UnsupportedInputExtension(String),
    #[error("unsupported output extension: {0}")]
    UnsupportedOutputExtension(String),
    #[error("gzip output requested, but output path is not gzip-compressed: {0}")]
    OutputPathMustBeGzip(PathBuf),
    #[error("compressed input {input} requires a gzip output path; got {output}")]
    CompressedInputRequiresGzipOutput { input: PathBuf, output: PathBuf },
    #[error("compressed stdin requires a gzip output path; got {0}")]
    CompressedStdinRequiresGzipOutput(PathBuf),
    #[error("invalid isoforms row in {path}:{line}: {message}")]
    InvalidIsoformsLine {
        path: PathBuf,
        line: usize,
        message: String,
    },
    #[error("duplicate transcript mapping in {path}:{line}: {transcript}")]
    DuplicateIsoformMapping {
        path: PathBuf,
        line: usize,
        transcript: String,
    },
    #[error("failed to build thread pool: {0}")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError),
    #[error("reader error: {0}")]
    Reader(#[from] genepred::reader::ReaderError),
    #[error("logger error: {0}")]
    Logger(#[from] log::SetLoggerError),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}
