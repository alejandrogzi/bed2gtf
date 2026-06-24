# Changelog

## v2.1.0 - 2026-06-24

### Fixed

- Updated `genepred` from `0.0.11` to `0.0.15`, pulling in the GFF fix that
  prevents unmapped transcripts from reusing their own ID as the parent gene ID.
- Fixed GFF output for unmapped transcripts so `mRNA` records no longer form
  self-parent loops; generated gene IDs now use the distinct `gene-<tx>` form.
- Collapsed multiple isoforms from the same gene into a single emitted `gene`
  feature instead of writing one gene row per transcript.
- Recomputed emitted gene coordinates from the union of all isoform spans, so
  the gene feature covers the full gene rather than only the leader transcript.
- Preserved deterministic input order while parsing and rendering records in
  parallel chunks.
- Added regression coverage for GTF/GFF gene-line de-duplication, cross-window
  isoform aggregation, span correction, and GFF self-parent prevention.

### Changed

- Bumped the crate version from `2.0.0` to `2.1.0`.
- Refactored conversion rendering into a two-pass process: first gather records
  and gene spans, then render ordered output windows.

## v2.0.0 - 2026-03-19

### Changed

- Reworked the converter around `genepred`.
- Added support for all BED input layouts handled by the current CLI.
- Added GFF/GFF3 output support alongside GTF.
- Added compressed input and output handling through the new conversion path.
- Updated the README, project logo, and Nextflow module assets.
- Added container and publish workflow support.
- Removed the older implementation in favor of the new v2 architecture.

## v1.9.3 - 2024-11-20

### Fixed

- Fixed issue #11.
- Expanded gzip reader support.
- Reported defective runtime lines through `log::error`.

## v1.9.2 - 2024-04-12

### Added

- Added `--no-gene` support.
- Added maximum memory usage reporting.

### Changed

- Refactored `BedRecord` handling.
- Moved toward owned CLI module state.

## v1.9.1 - 2024-02-20

### Fixed

- Disabled invalid UTR output paths.
- Added coordinate guards to avoid invalid feature coordinates.
- Fixed automatic version update handling.

## v1.9.0 - 2024-01-04

### Fixed

- Fixed exon numbering, including reverse-strand indexing.
- Fixed `--gz` flag handling.

### Changed

- Reworked buffered output writing through boxed writers.

## v1.8.0 - 2023-11-27

### Changed

- Rebuilt the `BedRecord` model.
- Added parallel conversion handlers.
- Fixed gene coordinate handling in the parallel path.
- Added a line-builder module.
- Improved isoform parsing with partial `memchr` use.

### Fixed

- Added `PartialEq` support for `BedRecord`.
- Fixed typos and documentation drift around the v1.8 changes.

## v1.7.0 - 2023-10-17

### Added

- Added parsing error handling.
- Added a pre-sorting step and early sorting support.

### Changed

- Updated documentation for the sorting workflow.

## v1.6.0 - 2023-10-02

### Fixed

- Reported missing isoforms instead of unwrapping and aborting.

## v1.5.0 - 2023-09-13

### Changed

- Updated output feature naming, including `five_prime_utr`.
- Added logging improvements and new success messaging.
- Added basic test modules.
- Updated clap metadata.

## v1.4.0 - 2023-08-19

### Fixed

- Fixed missing spaces in `exon_number` and ID attribute output.

## v1.3.0 - 2023-08-18

### Fixed

- Fixed trailing spaces in transcript and UTR lines.

## v1.2.0 - 2023-08-15

### Changed

- Reimplemented the core structs and conversion model.
- Removed the old implementation.
- Added project license metadata.

## Initial history - 2023-08-02

### Added

- Bootstrapped the `bed2gtf` project.
- Added initial Cargo metadata, README content, `.gitignore`, and isoform
  documentation updates.
