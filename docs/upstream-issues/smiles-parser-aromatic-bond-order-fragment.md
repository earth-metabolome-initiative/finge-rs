# smiles-parser: retained kekule bond order on aromatic bonds leaks into fragment canonicalization

Status: FIXED upstream on 2026-06-08. The report and tests below are kept for
reference. After the fix, finge-rs MAP4 is representation-independent: the
`map4_normalized_path_matches_raw` test confirms the raw input-order graph and
the RDKit-normalized graph produce byte-identical shingle sets across the whole
corpus.

Discovered while validating finge-rs MAP4. This was the remaining blocker to a
representation-independent (truly canonical) MAP4: the rooted circular-substructure
label of an atom depended on which equivalent representation of the parent
molecule it was computed from.

## Summary

`canonicalize()` and `render_rooted()` treat the kekule single/double order
stored on an aromatic bond as significant. For a closed aromatic ring this is
harmless (the whole-molecule canonical form is the same either way). But when a
radius-`r` environment is carved out of an aromatic ring, the fragment is a
broken, non-closable aromatic system, and the retained kekule orders leak into
its canonical rendering. Two representations of the same molecule that differ
only in retained kekule orders then produce different fragment labels.

This makes `rooted_environment_smiles` (and therefore MAP4 shingles)
representation-dependent: the raw-parsed graph and the kekulize-then-reperceive
graph of the same molecule give different labels for the same atom, even though
the two graphs are canonically identical as whole molecules.

## Reproducer

p-xylene, `Cc1ccc(C)cc1`, radius-2 environment of atom 6 (a ring CH carbon).

```rust
use smiles_parser::smiles::{AromaticityPolicy, Smiles};

let raw: Smiles = "Cc1ccc(C)cc1".parse().unwrap();
// An equivalent representation: kekulize, then re-perceive RDKit-default
// aromaticity (exactly what finge-rs does to normalize a molecule).
let norm = raw
    .kekulize_standalone()
    .unwrap()
    .perceive_aromaticity_for(AromaticityPolicy::RdkitDefault)
    .unwrap()
    .into_aromaticized();

// Whole-molecule canonical forms are identical:
assert_eq!(raw.canonicalize().render(), norm.canonicalize().render()); // "Cc1ccc(C)cc1"

// But the radius-2 environment label of atom 6 differs:
assert_eq!(raw.rooted_environment_smiles(6, 2, false).unwrap(),  "c(cc)c(C)c");
assert_eq!(norm.rooted_environment_smiles(6, 2, false).unwrap(), "c(c(C)c)cc"); // diverges
```

The two fragments are identical in atoms, edges, and local numbering. They differ
only in the bond order retained on the aromatic ring bonds:

```text
raw  fragment bonds (local):  (0,1) Single arom, (0,3) Single arom, (1,2) Single arom, (2,5) Single arom, (0,4) Single
norm fragment bonds (local):  (0,1) Single arom, (0,3) Double arom, (1,2) Double arom, (2,5) Single arom, (0,4) Single

raw  fragment canonicalize().render() = "Cc(c)ccc"   rooted_at(6) = "c(cc)c(C)c"
norm fragment canonicalize().render() = "Cc(ccc)c"   rooted_at(6) = "c(c(C)c)cc"
```

The raw parse stores aromatic bonds as `Single` + aromatic flag (no kekule
orders). The kekulize-then-reperceive round-trip keeps the kekule `Double`/`Single`
orders alongside the aromatic flag. `canonicalize()` lets those orders change the
output for the broken-ring fragment.

RDKit writes such broken-ring fragments with lowercase aromatic atoms and no
kekule double bonds, i.e. it matches the raw (`c(cc)c(C)c`) result. The norm
result is the wrong one.

## Root cause

The single/double order on a bond that is flagged aromatic is non-semantic: the
aromatic flag already determines the bond. Canonicalization should therefore
ignore (or normalize) the kekule order of aromatic bonds, so that two graphs
which are identical except for kekule assignment on aromatic bonds canonicalize
identically. Today they do not, and the difference becomes observable as soon as
the aromatic system is opened (as in a MAP4 environment fragment).

It is specifically not a numbering or branch-order problem: whole-molecule
`render_rooted(i)` is identical for raw and norm at every atom, the environment
bond sets are identical, and the fragment local numbering is identical. The only
differing input is the retained aromatic-bond order.

## Expectation

Representation invariance. For any molecule `M` and any representation-preserving
transform `T` (parse, `canonicalize`, kekulize + reperceive round-trip), and for
every atom `i` and radius `r`:

```text
M.rooted_environment_smiles(i, r, iso) == T(M).rooted_environment_smiles(i, r, iso)
```

Equivalently, at the canonicalization layer: two graphs that differ only in the
kekule single/double order assigned to aromatic bonds must produce equal
`canonicalize().render()` output, for fragments as well as whole molecules.

## Suggested fix direction

The kekule order on an aromatic bond must not influence canonical labeling or
rendering. Either:

1. Normalize at the canonicalization layer: when computing canonical invariants
   and emitting SMILES, treat an aromatic bond's order as a single canonical
   value (the aromatic flag), ignoring any retained `Double`/`Single`. This is
   the general fix and also covers `render_rooted` (which builds on the canonical
   ordering).
2. Or normalize at fragment construction: when `fragment_from_bonds` /
   `fragment_from_atoms` carve a fragment, reset the order of every aromatic bond
   to the canonical single+aromatic form (or re-kekulize the fragment to a
   canonical kekule structure). This is the narrower fix that covers the MAP4
   path specifically.

Option 1 is preferable because it makes canonicalization sound with respect to
kekule representation everywhere, not only on the fragment path.

## How to test

Representation invariance over a corpus of aromatic molecules:

```rust
use smiles_parser::smiles::{AromaticityPolicy, Smiles};

fn roundtrip(m: &Smiles) -> Smiles {
    m.kekulize_standalone()
        .unwrap()
        .perceive_aromaticity_for(AromaticityPolicy::RdkitDefault)
        .unwrap()
        .into_aromaticized()
}

for smiles in ["Cc1ccc(C)cc1", "Cc1ccccc1", "c1ccc2ccccc2c1", "O=C(O)c1ccccc1"] {
    let raw: Smiles = smiles.parse().unwrap();
    let norm = roundtrip(&raw);
    assert_eq!(raw.nodes().len(), norm.nodes().len());
    for i in 0..raw.nodes().len() {
        for r in 1..=3 {
            assert_eq!(
                raw.rooted_environment_smiles(i, r, false),
                norm.rooted_environment_smiles(i, r, false),
                "{smiles} atom {i} radius {r} is representation-dependent",
            );
        }
    }
}
```

Direct canonicalization invariance (the underlying property), comparing a raw
aromatic parse against its kekulized-then-reperceived form on a broken-ring
fragment:

```rust
let raw: Smiles = "Cc1ccc(C)cc1".parse().unwrap();
let norm = roundtrip(&raw);
let raw_frag = raw.atom_environment(6, 2).unwrap().to_fragment().unwrap();
let norm_frag = norm.atom_environment(6, 2).unwrap().to_fragment().unwrap();
assert_eq!(
    raw_frag.smiles().canonicalize().render(),
    norm_frag.smiles().canonicalize().render(),
);
```

Exact-value pin once fixed (both must equal the RDKit-style result):

```rust
assert_eq!(raw.rooted_environment_smiles(6, 2, false).unwrap(), "c(cc)c(C)c");
assert_eq!(norm.rooted_environment_smiles(6, 2, false).unwrap(), "c(cc)c(C)c");
```

## Impact on finge-rs MAP4

finge-rs currently sidesteps this by computing MAP4 on the raw input-order parsed
graph, which matches RDKit numbering and aromatic-bond representation, giving
99.6% exact pairwise-Tanimoto parity with RDKit MAP4. That path is validated by
`map4_shingle_tanimoto_parity` and `map4_environment_partition_parity`. The known
gap is tracked by `map4_normalized_path_high_self_similarity` (raw vs normalized
mean Jaccard 0.83). Once aromatic-bond order is normalized in canonicalization,
the raw and normalized paths converge, the normalized (RDKit-aromaticity) graph
becomes the production basis, and finge-rs MAP4 becomes representation-independent
and strictly more canonical than RDKit (which has its own order-dependent rooted
SMILES).
