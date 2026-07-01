# smiles-parser: `non_isomeric()` and fragment labels retain redundant brackets

Status: bug report with proposed fix, 2026-06-07. Discovered while validating
finge-rs MAP4 environment-partition parity against an RDKit oracle. This is the
"Class 2" divergence from that work: finge-rs over-brackets organic-subset atoms
in a way that breaks its own canonical invariant.

## Summary

`Smiles::non_isomeric()` clears the isotope and stereo of each atom but does not
re-evaluate whether the atom still needs brackets. The bracket-to-organic-subset
collapse is performed only inside `canonicalize()`, so any render path that does
not go through a fresh canonicalization keeps stale brackets. As a result:

```text
"C[13CH3]".canonicalize().non_isomeric().render()  ==  "C[CH3]"   (wrong)
"C[13CH3]".non_isomeric().canonicalize().render()  ==  "CC"       (right)
```

The output depends on the order of `non_isomeric()` vs `canonicalize()`, which
means `non_isomeric()` is not idempotent with respect to canonical form. Its own
doctest (src/smiles/mod.rs:1224) asserts
`labeled.non_isomeric().render() == plain.non_isomeric().render()`, which is
exactly the property this bug violates for isotope- or bracket-decorated atoms.

The same gap surfaces in the MAP4 fragment label path
(`Fragment::render_rooted` -> `non_isomeric().render_rooted(local)`): fragment
atoms inherit the parent's bracket form and are never re-collapsed, so a bracketed
aromatic nitrogen renders as `[n]` where it should be `n`.

## Reproducers (finge-rs vs RDKit 2026.03.2)

| input | finge `render()` | finge `canonicalize()` | finge `canon().non_isomeric()` | RDKit non-isomeric | correct? |
|---|---|---|---|---|---|
| `C[OH]`    | `C[OH]`    | `CO`        | `CO`      | `CO` | render not normalized, canon ok |
| `C[13CH3]` | `C[13CH3]` | `C[13CH3]`  | `C[CH3]`  | `CC` | **bug** |
| `[12CH4]`  | `[CH4]`    | `[12CH4]`   | `[CH4]`   | `C`  | **bug** |
| `C[NH3+]`  | `C[NH3+]`  | `C[NH3+]`   | `C[NH3+]` | `C[NH3+]` | ok (charge keeps bracket) |
| `c1cc[nH]c1` | `c1cc[nH]c1` | `c1cc[nH]c1` | `c1cc[nH]c1` | `c1cc[nH]c1` | ok (aromatic N-H) |
| `[2H]C`    | `[2H]C`    | `[2H]C`     | `[H]C`    | `[H]C` | ok (explicit H atom) |

MAP4 fragment labels (radius-1 rooted, non-isomeric), from
`[CH3][Ga]([CH3])[n]1cccn1`:

| atom | finge label | RDKit label |
|---|---|---|
| aromatic c next to ring N | `c(c)[n]` | `c(c)n` |
| ring N | `n(c)[n]` | `n(c)n` |

The symmetric case is the clearest invariant break: in `C[13CH3]` the two carbons
form one symmetric `C-C` fragment, yet finge-rs labels the plain carbon's
environment `CC` and the (isotope-stripped) carbon's environment `C[CH3]`. The
same fragment must get the same label.

## Root cause

Bracket-vs-organic-subset is **stored atom state** (`AtomSyntax::Bracket` vs the
organic-subset form built by `Atom::new_organic_subset`). `write_smiles`
(src/atom/mod.rs:555) emits whatever form the atom currently has. The collapse
from a bracket atom to the bare organic-subset form is done by
`maybe_collapse_atom_to_organic_subset` (src/smiles/canonicalization.rs:529),
which is called only during `canonicalize()`.

- `Atom::non_isomeric()` (src/atom/mod.rs:495) sets `isotope_mass_number = None`
  and `chirality = None` and returns. It cannot collapse, because the collapse
  needs molecule context (the atom's bonds, to compute the default implicit-H
  count) which an `Atom` does not have on its own.
- `Smiles::non_isomeric()` (src/smiles/mod.rs:1228) maps `atom.non_isomeric()`
  over the atoms and rebuilds the graph, but never runs the collapse pass even
  though it now has the molecule context to do so.
- `render()` / `render_rooted()` echo the stored form faithfully (by design:
  render is faithful, canonicalize normalizes). The MAP4 label path relies on
  `non_isomeric()` for normalization, and that step does not normalize brackets.

So once an atom is bracketed (because of an isotope in the input, or because the
parent forced a bracket that the fragment no longer needs), clearing the isotope
or carving a fragment leaves the bracket in place.

## Proposed fix

Make `Smiles::non_isomeric()` run the existing organic-subset collapse over every
atom after clearing isotope and stereo, reusing the canonicalization rule rather
than duplicating it. In `Smiles::non_isomeric()`, after the non-isomeric graph is
built:

```rust
// after building the non-isomeric `Self`, collapse atoms that no longer need
// brackets now that isotopes (and any other bracket-only reason) are gone.
let collapsed_atoms: Vec<_> = (0..result.atom_nodes.len())
    .map(|node_id| {
        maybe_collapse_atom_to_organic_subset(&result, node_id, result.atom_nodes[node_id])
    })
    .collect();
// rebuild `result` with `collapsed_atoms` (same bonds / sidecars)
```

`maybe_collapse_atom_to_organic_subset` already encodes the correct eligibility
rule, so the fix is to apply it, not to reinvent it. An atom collapses to the
bare symbol only when ALL of these hold (verbatim from the existing function):

- syntax is `Bracket`,
- `isotope_mass_number` is `None`,
- formal charge is `0`,
- atom-map class is `0`,
- chirality is `None`,
- not an aromatic element that cannot be written unbracketed
  (`can_write_unbracketed_aromatic`),
- symbol is in the unbracketed-valid set (`B C N O P S F Cl Br I` or wildcard),
- the implicit-H count it would get if written unbracketed equals its current
  hydrogen count.

Because this rule is exact, the controls in the table above stay bracketed:
`C[NH3+]` (nonzero charge), `c1cc[nH]c1` (aromatic N with an H count that differs
from the unbracketed default), `[H]C` (the explicit hydrogen is its own atom),
metals such as `[Ga]` (not in the organic subset), and anything with an atom map
or radical.

Since the MAP4 label path is `Fragment::render_rooted` ->
`non_isomeric().render_rooted(local)`, fixing `non_isomeric()` fixes both the
`C[13CH3]` case and the fragment `[n]` case in one place, and it restores
`non_isomeric()`'s own idempotence with respect to canonical form.

### Optional secondary normalization

`Fragment::render_rooted(parent, isomeric = true)` and a plain `render()` of a
non-canonicalized graph still echo stored brackets. That is consistent with the
"render is faithful, canonicalize normalizes" design, so it can be left as is.
If you would rather have fragment labels normalized regardless of the isomeric
flag, run the same collapse in `Fragment::render_rooted` before emitting, or
expose a small `Smiles::with_unbracketed_collapse()` helper that both
`non_isomeric()` and the fragment path call.

## Test cases to add

```rust
// non_isomeric() must be order-independent with respect to canonicalization
for s in ["C[13CH3]", "[12CH4]", "C[OH]"] {
    let m: Smiles = s.parse().unwrap();
    assert_eq!(
        m.canonicalize().non_isomeric().render(),
        m.non_isomeric().canonicalize().render(),
    );
}
assert_eq!("C[13CH3]".parse::<Smiles>().unwrap().non_isomeric().render(), "CC");
assert_eq!("[12CH4]".parse::<Smiles>().unwrap().non_isomeric().render(), "C");

// controls: these MUST keep their brackets
for s in ["C[NH3+]", "c1cc[nH]c1", "[H]C"] {
    let m: Smiles = s.parse().unwrap();
    assert!(m.non_isomeric().render().contains('['));
}

// MAP4 fragment label: bracket-forced aromatic N collapses in the fragment
let m: Smiles = "[CH3][Ga]([CH3])[n]1cccn1".parse().unwrap();
// the aromatic-carbon environment must read c(c)n, not c(c)[n]
assert!(
    m.rooted_environment_smiles(/* an aromatic c */, 1, false)
        .unwrap()
        .find("[n]")
        .is_none()
);
```

## Impact on finge-rs MAP4

This is the only divergence class that is a real defect (the other class is
RDKit's own order-dependent rooted SMILES, where finge-rs is the more canonical
side and no change is wanted). Fixing it removes every "consistency" violation in
the environment-partition probe, which makes finge-rs's MAP4 environment
partition a clean refinement of RDKit's: same fragment always gets the same
label, and labels differ from RDKit only where RDKit over-splits isomorphic
environments. After the fix we can quantify the residual finge-vs-RDKit Tanimoto
deviation as the headline MAP4 correctness number.
