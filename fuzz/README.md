# Fuzz Targets

This crate uses `cargo-fuzz` with three targets:

- `ecfp`
- `atom_pair`
- `topological_torsion`

Each target:

- accepts a fuzzed `String`
- parses it as a `smiles_rs::smiles::Smiles`
- returns early if parsing fails
- computes the relevant fingerprint(s)
- asserts basic structural invariants on the result

## Checked Invariants

### `ecfp`

- bit fingerprints have strictly increasing active-bit iterators
- all active bits are in range and readable through `contains()`
- counted fingerprints have strictly increasing active-count iterators
- all active counts are in range and nonzero
- bit ECFP active bits match the nonzero bins of counted ECFP
- layered counted ECFP sums back to counted ECFP
- the layered radius-`0` section matches `CountEcfpFingerprint::new(0, ..)`
- radius-`2` ECFP/count fingerprints are monotone extensions of radius-`1`
- the same invariants hold after `with_explicit_hydrogens()`

### `atom_pair`

- default and `countSimulation=false` fingerprints do not panic
- all active bits are in range and readable through `contains()`
- the same invariants hold on raw `Smiles`, explicit-H `Smiles`, and
  `SmilesRdkitScratch`-prepared graphs

### `topological_torsion`

- default and non-default toggles do not panic:
  - `countSimulation=false`
  - `onlyShortestPaths=true`
  - `torsionAtomCount=5`
- all active bits are in range and readable through `contains()`
- `onlyShortestPaths=true` produces a subset of the default prepared-graph bits
- the same basic invariants hold on raw `Smiles`, explicit-H `Smiles`, and
  `SmilesRdkitScratch`-prepared graphs

## Running

Install `cargo-fuzz` if needed:

```bash
cargo install cargo-fuzz
```

Then run a target, for example:

```bash
cargo fuzz run ecfp
cargo fuzz run atom_pair
cargo fuzz run topological_torsion
```

## Seed Corpus

Every target has a seed corpus in `fuzz/seeds/<target>/`, one input per file, and the ClusterFuzzLite build fails for a target without one. The seeds start from real SMILES in `tests/fixtures`, are grown by fuzzing, and are then reduced by libFuzzer's `-set_cover_merge=1 -use_counters=0` to the smallest set that keeps the same edge coverage. The targets that take a plain SMILES string also carry a few dozen hand-picked molecules covering charges, isotopes, stereo, aromatic rings and past crash inputs. The `mutator_*` targets decode their input through `arbitrary` as a SMILES string and a `u64` RNG seed, so their seeds are not readable SMILES. Pass the seeds as a second corpus directory to start a local run from them:

```bash
mkdir -p fuzz/corpus/ecfp
cargo fuzz run ecfp fuzz/corpus/ecfp fuzz/seeds/ecfp
```
