# finge-rs

[![CI](https://github.com/earth-metabolome-initiative/finge-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/earth-metabolome-initiative/finge-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Rust 1.86+](https://img.shields.io/badge/rust-1.86%2B-orange.svg)](Cargo.toml)

finge-rs computes molecular fingerprints in `no_std` Rust (with `extern crate alloc`), built around a small set of traits for molecular atoms, bonds, and graphs so any backend can supply the chemistry. Where it claims RDKit parity it reproduces RDKit's bit output exactly, checked against a tracked corpus of RDKit-parseable SMILES from the scikit-fingerprints HIV set across ECFP radii `0` to `5` and AtomPair, RDK, and Topological Torsion sizes from `64` to `4096`.

It provides the usual connectivity fingerprints. Morgan/ECFP comes as a bit fingerprint (`EcfpFingerprint`), a folded-count fingerprint (`CountEcfpFingerprint`), and an exact-radius folded-count variant (`LayeredCountEcfpFingerprint`), and alongside them are 2D AtomPair (`AtomPairFingerprint`), RDKit topological paths (`RdkFingerprint`), Topological Torsion (`TopologicalTorsionFingerprint`), and the 167-bit MACCS keys (`MaccsFingerprint`) behind the `smarts-support` feature. Each folds into a `BitFingerprint` or `CountFingerprint`.

MAP4 (`Map4Fingerprint`, MinHashed atom pair) rounds out the set. It exposes its raw shingle set, folds into a bit or counted fingerprint like the others, and produces a MinHash signature for similarity work, validated for Tanimoto parity against the reference `map4` implementation. ECFP can be sketched to a MinHash the same way. Those signatures feed `LshIndex`, a banded locality-sensitive-hashing index for approximate nearest-neighbour search over large collections, which can hand a k-NN graph straight to the [`bhtsne`](https://crates.io/crates/bhtsne) t-SNE behind the `tsne` feature.

`SmilesRdkitScratch` turns a `smiles-parser` molecule into an RDKit-normalized graph that every fingerprint accepts. AtomPair and Topological Torsion also run on a raw `Smiles` when you do not need the normalization step.

## Usage

The example below prepares a few molecules, folds ECFP and MAP4 fingerprints, and builds a MAP4 MinHash index to retrieve a molecule's nearest neighbours.

```rust
use finge_rs::{EcfpFingerprint, Fingerprint, LshIndex, Map4Fingerprint};
use finge_rs::smiles_support::SmilesRdkitScratch;
use smiles_parser::smiles::Smiles;

let mut scratch = SmilesRdkitScratch::default();
let mut index = LshIndex::<u32, 512>::new(256);

for smiles in ["CCO", "OCC1OC(O)C(O)C(O)C1O", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
    let molecule: Smiles = smiles.parse().expect("example SMILES should parse");
    let graph = scratch.try_prepare(&molecule).expect("preparation should succeed");

    let _ecfp = EcfpFingerprint::default().compute(&graph);
    let _map4_bits = Map4Fingerprint::default().compute(&graph);
    index.insert(Map4Fingerprint::default().minhash::<_, u32, 512>(&graph));
}

let query: Smiles = "CCO".parse().expect("example SMILES should parse");
let graph = scratch.try_prepare(&query).expect("preparation should succeed");
let signature = Map4Fingerprint::default().minhash::<_, u32, 512>(&graph);
let neighbours = index.query(&signature, 3);
assert!(!neighbours.is_empty());
```
