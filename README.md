# finge-rs

[![CI](https://github.com/earth-metabolome-initiative/finge-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/earth-metabolome-initiative/finge-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Rust 1.86+](https://img.shields.io/badge/rust-1.86%2B-orange.svg)](Cargo.toml)

Trait-first molecular fingerprints for `no_std` Rust with `extern crate alloc`.

Right now the crate provides:

- generic traits for molecular atoms, bonds, and graphs
- `BitFingerprint`
- `CountFingerprint`
- bit-only, non-chiral Morgan/ECFP through `EcfpFingerprint`
- folded-count, non-chiral Morgan/ECFP through `CountEcfpFingerprint`
- exact-radius folded-count Morgan/ECFP through `LayeredCountEcfpFingerprint`
- bit-only, non-chiral 2D AtomPair through `AtomPairFingerprint`
- bit-only, non-chiral Topological Torsion through `TopologicalTorsionFingerprint`
- optional RDKit-normalized `smiles-parser` integration behind `smiles-support`

Current RDKit parity coverage:

- 1024 `smiles-parser`- and RDKit-parseable SMILES from the tracked scikit-fingerprints HIV fixture corpus
- bit, folded-count, and exact-radius folded-count ECFP radii `0` through `5`
- ECFP, AtomPair, and Topological Torsion bit sizes `64`, `128`, `256`, `512`, `1024`, `2048`, and `4096`

## Usage

Under `smiles-support`, `SmilesRdkitScratch` prepares a RDKit-normalized
`smiles-parser` graph that works with both `EcfpFingerprint` and
`AtomPairFingerprint`/`TopologicalTorsionFingerprint`. `AtomPairFingerprint`
and `TopologicalTorsionFingerprint` can also run directly on raw
`Smiles` when you do not need the normalization step.

```rust
# #[cfg(feature = "smiles-support")]
# fn main() {
use finge_rs::{
    AtomPairFingerprint, CountEcfpFingerprint, EcfpFingerprint, Fingerprint,
    LayeredCountEcfpFingerprint,
    TopologicalTorsionFingerprint,
    smiles_support::SmilesRdkitScratch,
};
use smiles_parser::smiles::Smiles;

let smiles: Smiles = "CCO".parse().expect("example SMILES should parse");
let mut scratch = SmilesRdkitScratch::default();
let graph = scratch.try_prepare(&smiles).expect("fingerprint preparation should succeed");
let atom_pair = AtomPairFingerprint::default().compute(&graph);
let ecfp = EcfpFingerprint::default().compute(&graph);
let counted_ecfp = CountEcfpFingerprint::default().compute(&graph);
let layered_counted_ecfp = LayeredCountEcfpFingerprint::default().compute(&graph);
let torsion = TopologicalTorsionFingerprint::default().compute(&graph);
# let _ = atom_pair;
# let _ = ecfp;
# let _ = counted_ecfp;
# let _ = layered_counted_ecfp;
# let _ = torsion;
# }
# #[cfg(not(feature = "smiles-support"))]
# fn main() {}
```
