# RDKit Comparison Benchmarks

The Rust Criterion benchmarks and the RDKit Python benchmark use the same
default corpus:

```text
tests/fixtures/pubchem_benchmark_10000_smiles.txt.gz
```

This file contains 10,000 unique SMILES selected from the random PubChem corpus
that was already filtered to parse in both `smiles-parser` and RDKit. For a
fresh comparison corpus, generate a random PubChem sample and then keep only
SMILES that parse in both toolchains before benchmarking.

To generate a fresh `smiles-parser`-parseable sample from PubChem:

```bash
cargo run --release --features datasets --example sample_random_parseable_corpus -- \
  --output tests/fixtures/pubchem_benchmark_10000_smiles.txt.gz \
  --limit 10000
```

Run the Rust benchmarks:

```bash
cargo bench --features smiles-support --bench ecfp
cargo bench --features smiles-support --bench atom_pair
cargo bench --features smiles-support --bench topological_torsion
cargo bench --features smiles-support --bench rdk
cargo bench --features smarts-support --bench maccs
```

Run the matching RDKit benchmark:

```bash
uv run --with rdkit python benches/rdkit_fingerprints.py \
  --input tests/fixtures/pubchem_benchmark_10000_smiles.txt.gz \
  --fingerprint all
```

To use a different corpus for the Rust benches:

```bash
FINGE_RS_BENCH_SMILES=/path/to/smiles.txt.gz cargo bench --features smiles-support --bench ecfp
```
