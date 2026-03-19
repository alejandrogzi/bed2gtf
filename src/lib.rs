// Copyright (c) 2025 Alejandro Gonzalez-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

pub mod app;
pub mod cli;
pub(crate) mod config;
pub(crate) mod convert;
pub(crate) mod detect;
pub(crate) mod error;
pub(crate) mod isoforms;
pub(crate) mod memory;

pub use app::run_cli;
pub use cli::{Args, BedType, OutputFormat};
pub use error::{Bed2GtfError, Result};
