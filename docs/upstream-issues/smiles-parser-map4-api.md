# smiles-parser: methods needed to support the MAP4 fingerprint

Status: design proposal, 2026-06-06. Author-facing request for new `smiles-parser`
public methods. finge-rs is implementing MAP4 (MinHashed Atom Pair, diameter 4)
and needs three capabilities that `smiles-parser` does not currently expose.

## Why finge-rs needs this

MAP4 represents a molecule as a set of "shingles", one per (atom pair, radius).
Each shingle is the string

```text
min(CS_r(j), CS_r(k)) | topological_distance(j, k) | max(CS_r(j), CS_r(k))
```

where `CS_r(i)` is the canonical SMILES of the circular substructure of radius
`r` around atom `i`, rooted at `i`, written non-isomeric. RDKit builds `CS_r(i)`
in three steps:

1. `FindAtomEnvironmentOfRadiusN(mol, r, i)` returns the bond set of the
   radius-`r` circular environment of atom `i` (the same environment ECFP uses).
2. `PathToSubmol(mol, env, atomMap=...)` materializes that bond set as a
   standalone fragment and reports where atom `i` landed in the fragment.
3. `MolToSmiles(submol, rootedAtAtom=i_in_submol, canonical=True, isomericSmiles=False)`
   renders the fragment canonically, forcing the traversal to start at the
   center atom.

The shingle strings are never parsed back. They are opaque labels fed into a
MinHash. The fingerprint's meaning is the Jaccard overlap of two molecules'
shingle sets, which depends only on the equivalence relation "do these two
environments get the same label", not on the exact spelling of the label.

### This is explicitly NOT a byte-parity request

`smiles-parser`'s canonicalization is sound but does not produce byte-identical
output to RDKit, and it does not need to. What MAP4 requires of `CS_r(i)` is:

- Deterministic: same fragment plus same center always yields the same string.
- Center-anchored: the label distinguishes environments that are isomorphic as
  graphs but differ in where the center sits.
- Isomorphism-faithful: two environments get the same label if and only if they
  are isomorphic as center-rooted, typed (element, charge, bond order, aromatic
  vs aliphatic) graphs.

If those three hold, finge-rs MAP4 and RDKit MAP4 produce identical
Jaccard/Tanimoto for every molecule pair (the MinHash signatures differ by a
relabeling of the shingle vocabulary, but neighbor lists and similarity matrices
are identical). Byte parity is neither required nor pursued. The RDKit strings
shown below are the semantic reference, not a literal target. `smiles-parser`'s
own canonical dialect is fine as long as the contract above holds.

Because the label is only an equivalence key, "rooted rendering" can be
satisfied two ways, and either is acceptable:

- Route A: force the SMILES traversal to start at the center (RDKit's
  `rootedAtAtom`).
- Route B: give the center atom a distinguishing canonical invariant (a marker
  color) and take the ordinary free canonical form of the colored fragment.
  Coloring forces the center to map to itself under any isomorphism, so it
  induces the same equivalence relation as rooting. This reuses the existing
  `canonicalize()` path and may be far less work.

## The capability gap

| # | Need | Today in smiles-parser |
|---|------|------------------------|
| 1 | Radius-N circular environment of an atom (bond set + atom set) | absent (ECFP-style traversal exists in finge-rs, but the empty-shell semantics below are bug-prone and better owned here) |
| 2 | Build a fragment from a bond set or atom set, with an atom map and fragment-local re-perception | `induced_subgraph(&[usize]) -> Self` exists but is private, takes only an atom set, returns no atom map, and inherits the parent implicit-H cache instead of re-perceiving the fragment |
| 3 | Canonical render anchored on a chosen atom, non-isomeric | `render()` and `canonicalize()` pick their own root and have no isomeric toggle exposed for this path |

## Bug-prone semantics that must be pinned down

These are where reference MAP4 implementations diverge, and where finge-rs must
match RDKit's environment definition exactly. Measured against RDKit 2026.03.2.

### Empty-shell rule (most important)

`FindAtomEnvironmentOfRadiusN(mol, r, i)` returns EMPTY (not a truncated
environment) when the outermost requested shell adds no new bond, i.e. when the
center's eccentricity is less than `r`. When non-empty, it returns all bonds
from shells `1..=r` cumulatively.

```text
ethanol  CCO   atom 0 (terminal C):  r1 -> "CC"      r2 -> "CCO"
ethanol  CCO   atom 1 (central C):   r1 -> "C(C)O"   r2 -> ''   (eccentricity 1 < 2)
acetic   CC(=O)O  atom 1 (carbonyl): r1 -> "C(C)(=O)O"  r2 -> ''
methane  C     atom 0:               r1 -> ''   (no bonds at all)
[NH4+]         atom 0:               r1 -> ''
```

An empty environment yields an empty label, and MAP4 still emits a shingle for
it (for example `|2|CCO`). Empty labels are part of correct MAP4 output and must
not be skipped. finge-rs will handle the shingle emission, but the
environment/render methods must signal "empty" distinctly (return `None` or an
empty environment) rather than erroring or returning the center atom alone.

### Cumulative content and ring-closure bonds

A non-empty environment includes every bond within `r` bonds of the center,
including ring-closure bonds between two atoms that are both at the frontier.
This is what makes a benzene-ring environment render as a ring rather than a
chain. The bond set, not an atom-induced subgraph, is the faithful primitive:
an atom-induced subgraph can pull in a peripheral ring bond that the radius-`r`
bond set legitimately excludes.

### Fragment re-perception

RDKit sanitizes the submol as a standalone molecule: aromaticity is re-perceived
on the fragment and implicit-H counts are recomputed for the fragment's valences.
A fragment cut from an aromatic ring is usually a broken ring whose aromaticity
differs from the parent. The current `induced_subgraph` copies the parent's
`implicit_hydrogen_cache`, which is the wrong behavior for this use case. The new
subgraph method must re-perceive aromaticity and recompute implicit H on the
fragment. This is the single largest risk to equivalence parity and should be a
first-class, tested behavior.

### Non-isomeric output

MAP4 uses `isomericSmiles=False`: tetrahedral and double-bond stereo are dropped
and isotope labels are dropped. The rooted render must support this mode so that,
for example, `[13CH3]` and `[12CH3]` environments collapse to the same label.

## Proposed methods

Signatures are sketches against the current public types (`Smiles<AtomPolicy>`,
`BondEdge`, `SmilesCanonicalLabeling`, `usize` node ids). Exact names are the
maintainer's call.

### Method 1: `atom_environment`

Returns an opaque, borrowed view of the radius-`radius` circular environment of
`center`, matching RDKit `FindAtomEnvironmentOfRadiusN`. No owned collection is
handed out: the view borrows the parent, exposes iterators, and carries the
consumers MAP4 actually needs (fragment construction and rooted rendering).
Because the empty-shell case returns `None`, a live `AtomEnvironment` is always
non-empty and always knows its own center.

```rust
pub fn atom_environment(
    &self,
    center: usize,
    radius: usize,
) -> Option<AtomEnvironment<'_, AtomPolicy>>;

pub struct AtomEnvironment<'mol, AtomPolicy = ConcreteAtoms> { /* private */ }

impl<'mol, AtomPolicy: SmilesAtomPolicy> AtomEnvironment<'mol, AtomPolicy> {
    pub fn center(&self) -> usize;
    pub fn atom_count(&self) -> usize;
    pub fn bond_count(&self) -> usize;
    pub fn contains_atom(&self, atom_id: usize) -> bool;
    pub fn atoms(&self) -> impl Iterator<Item = usize> + '_;
    pub fn bonds(&self) -> impl Iterator<Item = BondEdge> + '_;

    /// Materialize as a standalone fragment (Method 2).
    pub fn to_fragment(&self) -> Result<Fragment<AtomPolicy>, SubgraphError>;
    /// One-shot MAP4 label: fragment, then render rooted at this environment's
    /// center. `isomeric == false` for MAP4.
    pub fn rooted_smiles(&self, isomeric: bool) -> Result<String, SubgraphError>;
}
```

Intended usage and expected output:

```rust
let mol: Smiles = "CCO".parse()?;
assert!(mol.atom_environment(1, 2).is_none());     // central C, eccentricity 1 < 2 -> empty
let env = mol.atom_environment(0, 2).unwrap();     // terminal C
assert_eq!(env.center(), 0);
assert_eq!(env.atom_count(), 3);                   // atoms 0, 1, 2
assert_eq!(env.bond_count(), 2);
let label = env.rooted_smiles(false)?;             // RDKit reference: "CCO"
```

Returns `None` for the empty-shell case and for a `center` with no incident
bonds. Errors:

- `center >= atom_count` -> this is a programming error. Prefer a panic with a
  clear message (consistent with other index-based accessors), or an
  `Err(EnvironmentError::AtomOutOfRange)` if the maintainer prefers fallible.
- `radius == 0` -> always `None` (an environment of radius 0 has no bonds). MAP4
  never requests radius 0, but the method should define it.

### Method 2: `Fragment`

A `Fragment` is a standalone molecule carved out of a parent, with the atom
correspondence kept private behind two lookups so callers never touch a raw map.
It is produced by `AtomEnvironment::to_fragment` (the MAP4 path) or directly from
a bond or atom set for other callers. Inputs are iterators, not slices, so no
`Vec` is required on either side of the call.

```rust
/// Build a fragment from a bond set (RDKit `PathToSubmol`), re-perceiving
/// aromaticity and recomputing implicit H on the fragment.
pub fn fragment_from_bonds(
    &self,
    bonds: impl IntoIterator<Item = BondEdge>,
) -> Result<Fragment<AtomPolicy>, SubgraphError>;

/// Atom-set variant (induced subgraph): includes every bond among the atoms.
pub fn fragment_from_atoms(
    &self,
    atoms: impl IntoIterator<Item = usize>,
) -> Result<Fragment<AtomPolicy>, SubgraphError>;

pub struct Fragment<AtomPolicy = ConcreteAtoms> { /* private */ }

impl<AtomPolicy: SmilesAtomPolicy> Fragment<AtomPolicy> {
    pub fn smiles(&self) -> &Smiles<AtomPolicy>;
    pub fn into_smiles(self) -> Smiles<AtomPolicy>;
    pub fn atom_count(&self) -> usize;
    /// Fragment-local id of a parent atom, if it is in the fragment.
    pub fn local_id(&self, parent_atom: usize) -> Option<usize>;
    /// Parent id of a fragment atom.
    pub fn parent_id(&self, local_atom: usize) -> usize;
    /// Canonical render rooted at a parent atom (Method 3), resolving the map
    /// internally. Errors if the parent atom is not in the fragment.
    pub fn render_rooted(&self, parent_atom: usize, isomeric: bool)
        -> Result<String, RootError>;
}
```

Intended usage and expected output:

```rust
let mol: Smiles = "CCO".parse()?;
let frag = mol.atom_environment(0, 2).unwrap().to_fragment()?;
assert_eq!(frag.atom_count(), 3);
let label = frag.render_rooted(0, false)?;   // parent coords; RDKit reference: "CCO"
```

Behavioral requirements:

- Re-perceive aromaticity on the fragment (do NOT inherit the parent cache).
- Recompute implicit-H counts for the fragment's valences.
- Disconnected atom/bond sets are allowed (canonicalization already handles
  multiple components). MAP4 environments are connected by construction.

Errors (`SubgraphError`):

- `AtomOutOfRange` / `BondReferencesUnknownAtom`: an id in the input is not a
  valid atom of `self`.
- `PerceptionFailed`: aromaticity perception or kekulization of the fragment
  failed. This can happen for pathological broken-ring fragments. MAP4 should be
  able to treat this as "no usable label" (empty), so a typed error that the
  caller can map to an empty label is preferable to a panic.

### Method 3: `render_rooted` and `canonical_labeling_rooted`

```rust
/// Canonical SMILES of `self`, traversal anchored at `root`.
/// `isomeric == false` drops stereo and isotopes (MAP4 uses false).
pub fn render_rooted(&self, root: usize, isomeric: bool) -> Result<String, RootError>;

/// The canonical labeling with `root` forced to ordinal 0, for callers that
/// want the ordering rather than a string.
pub fn canonical_labeling_rooted(&self, root: usize)
    -> Result<SmilesCanonicalLabeling, RootError>;
```

Intended usage and expected output (RDKit strings shown as the semantic
reference, smiles-parser's canonical spelling may differ but must be stable and
isomorphism-faithful):

```rust
let mol: Smiles = "CCO".parse()?;
// Cleanest path: ask the environment directly (composes Methods 1 to 3).
let label = mol.atom_environment(1, 1).unwrap().rooted_smiles(false)?;
// central C, r1; RDKit reference: "C(C)O"

let mol: Smiles = "Cc1ccccc1".parse()?;               // toluene
// atom 1 (ring C bearing the methyl), r2, RDKit reference: "c(C)(cc)cc"
// atom 0 (methyl C), r1, RDKit reference: "Cc"
```

If Route B (center coloring) is chosen instead of true rooting, this method can
be implemented as "color `root`, canonicalize, render", and the returned string
will differ from the RDKit reference while satisfying the same equivalence
contract. That is acceptable.

Errors (`RootError`):

- `root >= atom_count`: index error, panic-or-`Err` per the maintainer's
  convention for the rest of the crate.

### Method 4 (optional convenience): `rooted_environment_smiles`

```rust
/// One-shot MAP4 substructure label: environment -> fragment -> rooted render.
/// Returns `None` for the empty-shell case (RDKit's empty-string label).
pub fn rooted_environment_smiles(
    &self,
    center: usize,
    radius: usize,
    isomeric: bool,
) -> Option<String>;
```

This is the exact primitive MAP4 calls once per (atom, radius). It composes
methods 1 to 3 and folds both "empty environment" and "perception failed" into
`None`, which the caller renders as the empty label. If the maintainer would
rather keep `smiles-parser` low-level, finge-rs can compose 1 to 3 itself and
this method can be skipped.

## End-to-end composition (what finge-rs will do)

```rust
// CS_r(i) for every atom and radius in 1..=2, matching map4._get_atom_envs.
fn cs(mol: &Smiles, center: usize, radius: usize) -> String {
    mol.rooted_environment_smiles(center, radius, /* isomeric = */ false)
        .unwrap_or_default()      // empty-shell -> "" label
}
```

Distances come from finge-rs's own BFS (already used for AtomPair, including the
RDKit `1e8` disconnected sentinel), so no distance API is needed from
`smiles-parser`.

## Open questions for the maintainer

1. Route A (true rooted render) or Route B (center coloring plus existing
   `canonicalize`)? Route B looks cheaper and is equivalence-equivalent. Route A
   is closer to RDKit and generally useful beyond MAP4.
2. Index errors: panic (like other `*_by_id` accessors) or fallible `Result`?
   This proposal sketches fallible for the subgraph and render paths because
   fragment perception can genuinely fail at runtime, but matches your house
   style otherwise.
3. Is fragment-local re-perception acceptable to add to a public subgraph
   constructor, or should it be a separate `sanitized_subgraph` to avoid
   surprising existing callers of the private `induced_subgraph`?
4. Does the existing `render()` already have an isomeric/non-isomeric switch we
   can route through, or does that need to be added alongside `render_rooted`?
