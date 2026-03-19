// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::Parser;
use std::process::ExitCode;

use bed2gtf::{Args, Bed2GtfError};

fn main() -> ExitCode {
    match bed2gtf::run_cli(Args::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(Bed2GtfError::Logger(err)) => {
            eprintln!("ERROR: {err}");
            ExitCode::FAILURE
        }
        Err(_) => ExitCode::FAILURE,
    }
}
