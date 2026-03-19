// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use crate::config::{Config, InputSource, OutputTarget};
use crate::convert::run;
use crate::error::Result;
use crate::Args;
use log::{debug, error, info};

/// Runs the CLI application with the given arguments.
///
/// ```rust,ignore
/// # use bed2gtf::Args;
/// let args = Args::parse_from(["bed2gtf", "--to", "gtf"]);
/// bed2gtf::run_cli(args)?;
/// # Ok::<(), bed2gtf::Bed2GtfError>(())
/// ```
pub fn run_cli(args: Args) -> Result<()> {
    simple_logger::init_with_level(args.level)?;
    debug!("DEBUG: Parsed args: {args:?}");

    let config = match Config::from_args(args) {
        Ok(config) => config,
        Err(err) => {
            error!("ERROR: {err}");
            return Err(err);
        }
    };

    match &config.input {
        InputSource::Path(path) => {
            info!("INFO: Reading input from {}", path.path.display());
        }
        InputSource::Stdin => {
            info!("INFO: Reading input from stdin");
        }
    }

    match &config.output {
        OutputTarget::Path(path) => {
            info!("INFO: Writing output to {}", path.path.display());
        }
        OutputTarget::Stdout => {
            info!("INFO: Writing output to stdout");
        }
    }

    info!("INFO: Using {} threads", config.threads);
    debug!("DEBUG: Resolved config: {config:?}");

    let stats = match run(&config) {
        Ok(stats) => stats,
        Err(err) => {
            error!("ERROR: {err}");
            return Err(err);
        }
    };

    info!("INFO: Elapsed: {:.4} secs", stats.elapsed.as_secs_f32());
    info!("INFO: Memory: {:.2} MB", stats.mem_delta_mb);
    Ok(())
}
