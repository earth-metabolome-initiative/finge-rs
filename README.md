# finge-rs

[![CI](https://github.com/earth-metabolome-initiative/finge-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/earth-metabolome-initiative/finge-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Rust 1.86+](https://img.shields.io/badge/rust-1.86%2B-orange.svg)](Cargo.toml)

Trait-first molecular fingerprints for `no_std` Rust with `extern crate alloc`.

Right now the crate provides:

- generic traits for molecular atoms, bonds, and graphs
- `BitFingerprint`
- bit-only, non-chiral Morgan/ECFP through `EcfpFingerprint`
- optional `smiles-parser` integration behind `smiles-support`

Current RDKit parity coverage:

- 128 parser-compatible molecules from the `scikit-fingerprints` HIV corpus
- radii `0` through `5`
- bit sizes `64`, `128`, `256`, `512`, `1024`, `2048`, and `4096`

## Usage

For repeated `Smiles` fingerprinting, reuse `SmilesEcfpScratch` and call the
normal `Fingerprint::compute` API on the prepared graph.

```rust
# #[cfg(feature = "smiles-support")]
# fn main() {
use finge_rs::{EcfpFingerprint, Fingerprint, smiles_support::SmilesEcfpScratch};
use smiles_parser::smiles::Smiles;

let smiles: Smiles = "CCO".parse().expect("example SMILES should parse");
let mut scratch = SmilesEcfpScratch::default();
let graph = scratch.prepare(&smiles);
let fingerprint = EcfpFingerprint::default().compute(&graph);
# let _ = fingerprint;
# }
# #[cfg(not(feature = "smiles-support"))]
# fn main() {}
```
