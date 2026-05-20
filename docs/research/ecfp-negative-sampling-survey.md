# ECFP-Detectable Chemistry Violations as an Auxiliary Negative-Sampling Task: A Literature Survey

## Summary

The exact combination proposed for `finge-rs` — generating chemically invalid molecules whose violations are *detectable from ECFP atom/bond invariants*, computing ECFP, and training an auxiliary per-violation-class classifier head — does not appear to have been published. The closest prior art comes from four threads: (i) molecular contrastive-learning frameworks (MolCLR, iMolCLR, GraphCL) that perturb molecules as augmentations and now treat "faulty negatives" as a *problem to mitigate* rather than a signal to exploit; (ii) SMILES-validity research (UnCorrupt SMILES / Schoenmaker et al.; "Invalid SMILES are beneficial...") which is the only line of work that explicitly classifies parse/sanitization failures into a discrete taxonomy of error types; (iii) generative-model evaluation (MOSES, GuacaMol) which defines validity but treats it as a scalar metric, not a multi-class label; and (iv) reward-shaped RL such as TSSR which separately scores valence, aromaticity, and connectivity violations. Reymond's GDB enumerations provide the canonical "rules of chemical validity" that any taxonomy can be anchored against. The clear gap: nobody has trained on a fingerprint (let alone ECFP specifically) with an explicit multi-class violation-type auxiliary head; the validity-as-auxiliary-signal idea has been tested on GNNs and language models, never on the count/bit fingerprint that drug-discovery production pipelines actually use.

> **Scope note.** `finge-rs` computes ECFP directly from its own graph types (`MolecularGraph` and trait impls) in pure Rust. RDKit appears in this document only as a reference point — (a) as the bit-identity target for the ECFP hash output, and (b) as the toolkit used by *prior published work* whose mutation strategies we are reviewing. No part of the `finge-rs` runtime invokes RDKit or its sanitizer.

---

## 1. What has actually been tried?

### 1a. Hard-negative / faulty-negative molecules in contrastive and self-supervised settings

- **MolCLR** (Wang, Wang, Cao, Farimani, 2021/2022): canonical molecular contrastive-learning framework. Uses three graph augmentations — atom masking, bond deletion, subgraph removal — that *can* produce chemically nonsensical graphs as positive pairs of the same anchor. Negatives are other molecules in the batch (no explicit invalidity construction). arXiv:2102.10056; *Nat. Mach. Intell.*, https://www.nature.com/articles/s42256-022-00447-x.
- **iMolCLR** (Wang et al., 2022): the direct follow-up. Explicitly identifies "faulty negatives" — pairs that should be similar but are pushed apart — and *down-weights* them using fingerprint (Tanimoto on ECFP) similarity. Adds fragment-level contrast. Closest published use of ECFP-like fingerprints in a contrastive-negative role, but treats fingerprint similarity as a *correction* to the negative set, not as a label for explicit invalidity. arXiv:2202.09346; *JCIM* 2022.
- **GraphCL** (You et al., 2020): node dropping, edge perturbation, attribute masking, subgraph as graph augmentations; chemistry validity not enforced. arXiv:2010.13902.
- **3D-Mol** (Li et al., 2023): topological-fingerprint dissimilarity weights negatives in 3D contrastive learning. arXiv:2309.17366.
- **DIG-Mol** (2024): dual-interaction GNN with momentum distillation; same paradigm. arXiv:2405.02628.

None train an explicit "is this molecule valid?" head, none explicitly *construct* invalids as a separate class.

### 1b. Explicit "is this molecule chemically valid?" classifiers

- **UnCorrupt SMILES** (Schoenmaker, Liu, Olivecrona, van Westen, Bender, 2023): transformer that *repairs* invalid SMILES. Pipeline classifies parse errors into a discrete set of types and uses those for analysis. Code: https://github.com/LindeSchoenmaker/SMILES-corrector. *J. Cheminf.* 15:33, https://link.springer.com/article/10.1186/s13321-023-00696-x.
- **"Invalid SMILES are beneficial rather than detrimental to chemical language models"** (Skinnider, *Nat. Mach. Intell.* 2024): parsing errors sorted into **six error types**; Cohen's d computed between loss distributions of valid vs each invalid type. Cleanest published per-violation-type taxonomy located. Counter-intuitive headline: forcing models toward validity hurts learning of other chemical properties. https://www.nature.com/articles/s42256-024-00821-x.
- **TSSR** (arXiv:2601.04521): per-category reward — `valence`, `aromaticity`, `connectivity` scored separately. Essentially the multi-class-violation idea but used as RL reward shaping inside a SMILES generator, not as an auxiliary classification head on a fixed encoder.
- **SmiSelf** (arXiv:2509.23099): separates syntactic from semantic invalidity, repairs via SELFIES round-trip. No discriminator.
- **MolGAN** (Cao & Kipf, 2018, arXiv:1805.11973): GAN discriminator distinguishes generated vs real molecular graphs — implicitly binary validity signal tangled with distributional realism.

### 1c. Sanitization failures as supervised labels

Schoenmaker's corrector and the Skinnider paper are the only works using sanitization failures *as labels* (not as a filter). Most other works use sanitization only to **discard** invalid samples (MOSES, GuacaMol — §4).

### 1d. Counterfactual / perturbed molecules via attribute or topology edits

- **MMACE** (Wellawatte, Seshadri, White, 2022, *Chem. Sci.* https://doi.org/10.1039/D1SC05259D): minimum-edit counterfactuals that flip a classifier's decision.
- **MMGCF** motif-rebuild counterfactual GNN framework, https://www.scirp.org/journal/paperinformation?paperid=140395.
- **XPlore / INDUCE / LLM-GCE**: graph counterfactual explanation methods (arXiv:2603.04209, 2306.04835, 2410.15165). All emphasize *valid* counterfactuals — invalids treated as failures, not training signal.
- **MMP / `mmpdb`** (Dalke, Hert, Kramer, 2018, https://pubs.acs.org/doi/10.1021/acs.jcim.8b00173): formal database of single-transformation edits, all valid by construction.

---

## 2. Taxonomies of invalidity

No canonical citable taxonomy of "ECFP-detectable chemistry violations" exists. What is available, ranked:

1. **Skinnider 2024 (six error types)** — closest to what `finge-rs` needs. Categories derived from chemistry-toolkit exception messages (kekulization, aromaticity, valence, ring closure, parenthesis/syntax, atom type). *Nat. Mach. Intell.* 2024.
2. **TSSR (arXiv:2601.04521)** — three buckets used as separate reward terms: **valence**, **aromaticity**, **connectivity**.
3. **Syntactic vs. semantic** (SmiSelf, SELFormer arXiv:2304.04662) — coarse two-way split across SMILES-generation literature. Out of scope for ECFP because ECFP cannot see SMILES syntax.
4. **Classical chemistry textbook taxonomy** — hypervalent / hypovalent / radical / charge-imbalance. Not from an ML paper; see Cooper et al., *Chem. Sci.* 2015 https://pubs.rsc.org/en/content/articlehtml/2015/sc/c5sc02076j for a quantitative definition of hypervalency.
5. **GDB construction rules** (Reymond group): valence + chemical stability + synthetic feasibility. The implicit "rules a molecule must satisfy" enumerated for GDB-11/13/17 (https://pubs.acs.org/doi/10.1021/ci300415d) is the closest *positive* taxonomy of validity; negate each rule to derive a violation class.

**Recommendation for `finge-rs`:** combine Skinnider's six-bucket structure with TSSR's three reward terms and *intersect with ECFP-detectability*. The resulting taxonomy, aligned with the atom-invariant tuple `(Z, degree, nH, q, Δmass, ring)` and the bond invariant, gives roughly:

| # | Violation class | ECFP signal channel |
|---|---|---|
| 1 | Hypervalent atom (degree > maxValence(Z)) | `degree`, optionally `nH` |
| 2 | Impossible formal charge for element | `q` |
| 3 | Nonexistent / implausible isotope (`Δmass` not in known mass table) | `Δmass` |
| 4 | Aromatic flag on atom in non-cycle / non-planar context | `ring`, bond type = 12 |
| 5 | Bond-order / valence mismatch (Σ bond orders ≠ effective valence + nH + q) | bond invariant + `degree`/`nH` |
| 6 | Unknown atomic number (Z=0 or > 118) | `Z` |
| 7 | Topological pathology (atom flagged in-ring but degree-1; disconnected fragment marked as ring member) | `ring` × `degree` |
| 8 | Ring aromaticity inconsistency (some bonds aromatic, others not, in same SSSR ring) | bond-invariant heterogeneity inside a ring |

This is novel synthesis — no single paper publishes this matrix.

---

## 3. Generation strategies — perturb vs. sample-from-broken

Two paradigms exist:

- **Perturb a real molecule (preferred).** All counterfactual-explanation work (MMACE, MMGCF, XPlore), all MMP-style methods (`mmpdb`), all evolutionary-mutation work (Acosta Murillo et al. 2025, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12692386/ — benchmarking graph-GA, graph-generative, SmilesClickChem, SELFIES-token, SMILES-token mutators), and NoiseMol (https://www.sciencedirect.com/science/article/abs/pii/S1093326323000529). **Why preferred:** keeps the molecule in a recognisable region of chemical space and isolates *one* violation per sample, giving clean labels.
- **Sample from a broken-molecule distribution (rare).** GDB-style enumeration of valence-violating graphs appears in Reymond's work only as an anti-example (filtered out). MolGAN's failure modes implicitly produce broken samples but they are unlabeled.

**Implementation patterns observed in prior (Python) code:**

- Graph-level edits via `RDKit.Chem.RWMol` (`AddAtom`, `RemoveAtom`, `AddBond`, `RemoveBond`, `ReplaceAtom`). Standard pattern documented in RDKit *Getting Started* and used by Schoenmaker's corrector training-data generator. Cited here as **prior-art reference only** — not applicable to the `finge-rs` runtime.
- SMILES token-level edits (mask / swap / delete) — NoiseMol; the "Going beyond SMILES enumeration" 2025 paper (RSC *Digital Discovery*, https://pubs.rsc.org/en/content/articlehtml/2025/dd/d5dd00028a) compares atom-deletion / random-masking / bioisosteric / functional-group masking at varying perturbation probabilities.
- SELFIES-token edits — guaranteed *syntactically* valid; bypass parser-level corruption (Krenn et al. 2020, *Mach. Learn.: Sci. Tech.*, arXiv:1905.13741). **Useless for `finge-rs`** because SELFIES eliminates the kinds of violation we *want* to detect.

For ECFP-detectable violations, **direct edits to per-atom / per-bond invariant entries** are the right primitive — they map directly to ECFP atom-invariant entries. In `finge-rs` this means mutating the underlying `MolecularGraph` fields (`atomic_numbers`, `hydrogen_counts`, `formal_charges`, `in_ring`) and the `ValuedCSR2D` bond-type entries.

---

## 4. Benchmarks and datasets

There is **no standard labeled corpus of invalid molecules**. Standard *valid*-molecule benchmarks are used:

- **MOSES** (Polykovskiy et al., *Front. Pharmacol.* 2020, https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.565644/full, code https://github.com/molecularsets/moses) — ZINC Clean Leads subset (~4.6 M). Reports **Validity** (parser-based), **Uniqueness**, **Novelty**, **Filters** (PAINS + MCF). All scalar.
- **GuacaMol** (Brown, Fiscato, Segler, Vaucher, *JCIM* 2019, https://pubs.acs.org/doi/10.1021/acs.jcim.8b00839) — ChEMBL-based generative benchmark with validity, KL-divergence, FCD, goal-directed sub-benchmarks. Scalar validity.
- **GDB-11/13/17** (Reymond group) — canonical *positive* set of chemically-valid molecules of bounded size. Useful as clean-positive reservoir.
- **ZINC / ChEMBL** — standard positives used by everyone.

Closest to a labeled-invalid corpus is the **SMILES-corrector training data on Zenodo (record 7157412)** — synthetic invalid SMILES paired with valid originals, labelled by error category. Hand-built for that paper, not adopted as a community benchmark.

**Replicable construction recipe (Schoenmaker et al.):** sample valid SMILES → token-level corruption (insertion, deletion, substitution, ring-number perturbation, parenthesis-balance corruption) → run through parser → keep ones that fail → record exception class. For `finge-rs` replace token-corruption with graph-level invariant-tuple corruption.

---

## 5. Gaps — specifically the combination "ECFP + auxiliary validity head + per-violation-class taxonomy"

**Not done.** No paper located that:

1. Computes ECFP / Morgan / circular fingerprints on **invalid** molecules and uses them as inputs.
2. Trains a multi-class classification head where classes correspond to *types* of chemical violation.
3. Couples this with a main property/affinity head as an auxiliary loss.

Closest prior art, ranked by proximity:

- **Skinnider 2024** — six-class violation taxonomy, but applied to language-model loss analysis, not a discriminative head, on tokens not fingerprints.
- **TSSR** (arXiv:2601.04521) — per-category violation rewards (valence/aromaticity/connectivity), used in RL on SMILES strings inside a generator, not an auxiliary head on a fixed encoder.
- **iMolCLR** — uses ECFP similarity to weight negatives in contrastive learning, but its negatives are *valid* molecules.
- **GROVER** (Rong et al., NeurIPS 2020, arXiv:2007.02835) — multi-task SSL pretraining with atom/bond context prediction and motif prediction as auxiliary heads. Same *architectural* pattern (encoder + multiple auxiliary heads) but auxiliary tasks are properties of *valid* molecules.
- **Hu et al.** *Strategies for Pre-training GNNs* (ICLR 2020, arXiv:1905.12265) — Context Prediction + Attribute Masking + graph-level supervised auxiliary prediction. Foundational template for auxiliary heads, never with a *validity-violation* head.
- **SCAGE** (*Nat. Commun.* 2025, https://www.nature.com/articles/s41467-025-59634-0) — four-task multitask pretraining including **molecular fingerprint prediction** as one auxiliary task. Reverse direction from `finge-rs`: predicts fingerprints from an encoder rather than treating fingerprints as inputs.
- **MolGAN** discriminator — implicitly binary validity head, no per-class taxonomy, no ECFP.

The unique design space `finge-rs` would occupy: **(fingerprint input) × (multi-class violation head) × (per-invariant-channel taxonomy)**. None of the three dimensions individually is novel; the intersection appears to be open.

---

## 6. Implementation patterns

### What the literature actually does, with code-level specifics

Listed as **prior-art reference for how Python pipelines have done this**, not as a template for `finge-rs`:

- **`RDKit.Chem.RWMol` graph edits** (Schoenmaker, MMP tools, MMACE). Typical Python pattern: load via `Chem.MolFromSmiles`, wrap in `RWMol`, mutate atom/bond properties (`SetAtomicNum`, `SetFormalCharge`, `SetIsotope`, `SetNumExplicitHs`, `AddBond`, `RemoveBond`), then selectively disable sanitization steps via `Chem.SanitizeMol(..., sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)` so the invalid molecule survives. Record which sanitization op would have failed → that is the violation label. Documented at https://www.rdkit.org/docs/GettingStartedInPython.html.
- **SMILES token-level corruption** (UnCorrupt SMILES, NoiseMol). Regex-tokenize a SMILES, mask / swap / delete tokens with probability `p`. Not relevant to `finge-rs`: most token corruptions either fail parsing (no ECFP computable) or produce a different *valid* molecule.
- **SELFIES-token corruption.** Cannot produce ECFP-detectable invalidity (SELFIES is 100% valid). Krenn et al., arXiv:1905.13741.
- **Graph-augmentation as data-loader transforms** (MolCLR, GraphCL). PyTorch-Geometric `Data` objects with `atom_masking`, `bond_deletion`, `subgraph` ops.

### Code-level guidance for `finge-rs`

Because the hard constraint is *ECFP-detectability*, mutations should be specified in the invariant-tuple space, not the SMILES space. In `finge-rs` this means direct edits to the `MolecularGraph` fields. Each operator is a localised edit:

| Mutation operator | Atom-invariant channel hit | `finge-rs` data structure to edit |
|---|---|---|
| `swap_atomic_number(i, Z')` | `Z` | `atomic_numbers[i]` |
| `swap_isotope(i, A')` | `Δmass` | (isotope vector — to add) |
| `set_formal_charge(i, q')` | `q` | `formal_charges[i]` |
| `add_excess_bond(i, j)` | `degree`, bond type | `edges` (ValuedCSR2D) |
| `delete_required_h(i)` | `nH` | `hydrogen_counts[i]` |
| `mark_aromatic_nonring(i)` | `ring` × bond=12 | `in_ring[i]`, bond entries |
| `break_aromatic_ring(i, j)` | `ring` × bond invariant | bond entries |

This sidesteps the sanitization gymnastics that every Python pipeline fights with: there is no sanitizer in the `finge-rs` runtime to fight, and validity is whatever the operator-set defines.

---

## Gaps and opportunities — concrete suggestions for `finge-rs`

1. **Publish the taxonomy.** The 8-row matrix in §2 (each row aligned to one ECFP invariant channel) does not exist in the literature in this exact form. Even ignoring the negative-sampling application, the matrix itself is a citable contribution.

2. **Build the operators in Rust on `MolecularGraph` directly.** All prior work routes mutations through Python + RDKit's `RWMol` + selective sanitization toggles. `finge-rs` can implement them as pure graph-level mutators on the existing per-atom vectors and the `ValuedCSR2D` bond store — faster, reproducible, and free of any external chemistry-toolkit dependency on the negative-generation path.

3. **Couple operators to invariant channels by construction.** Each operator targets exactly one invariant-tuple entry. This makes per-class labels mechanical and the negative *guaranteed* to be ECFP-detectable. No published work has this property.

4. **Two corpora, not one.** Following MOSES/GuacaMol conventions, ship (a) a positive corpus drawn from a clean source (ZINC / ChEMBL / GDB-13 subset), and (b) a paired-negative corpus where each negative is a single-operator perturbation of a positive with a recorded class label. This is the MMP convention applied to invalidity.

5. **Auxiliary head design.** Closest published template is GROVER (arXiv:2007.02835) — encoder + multiple auxiliary classification heads with weighted losses. For `finge-rs`, the auxiliary head is **multi-label** (a molecule can hit multiple violation classes simultaneously, e.g. hypervalent + aromaticity mismatch) — published work has only treated multi-class single-label.

6. **Evaluate against the Skinnider 2024 finding.** That paper claims invalid SMILES *help* language models. For ECFP + downstream classifier the question is open: does an invalidity-aware representation generalise better on standard MoleculeNet / Therapeutics Data Commons tasks? Publishable comparison.

7. **Use iMolCLR as the main comparison baseline.** Only published method using ECFP-derived similarity to shape contrastive negatives. `finge-rs` would generalise that "ECFP knows what's similar" intuition to "ECFP knows what's *wrong*."

8. **Avoid the SELFIES trap.** SELFIES makes the negative-construction problem trivially impossible. Any user wanting to integrate `finge-rs` negatives with a SELFIES-based pipeline needs an explicit warning in the crate docs.

---

## Short bibliography

- Wang Y., Wang J., Cao Z., Barati Farimani A. *MolCLR.* *Nat. Mach. Intell.* 2022. arXiv:2102.10056.
- Wang Y., Magar R., Liang C., Barati Farimani A. *iMolCLR.* *JCIM* 2022. arXiv:2202.09346.
- Schoenmaker L., Béquignon O.J.M., Jespers W., van Westen G.J.P. *UnCorrupt SMILES.* *J. Cheminf.* 15:33, 2023. https://link.springer.com/article/10.1186/s13321-023-00696-x. Code: https://github.com/LindeSchoenmaker/SMILES-corrector.
- Skinnider M.A. *Invalid SMILES are beneficial rather than detrimental to chemical language models.* *Nat. Mach. Intell.* 2024. https://www.nature.com/articles/s42256-024-00821-x.
- TSSR authors. *Two-Stage Swap-Reward-Driven RL for Character-Level SMILES Generation.* arXiv:2601.04521.
- SmiSelf authors. *How to Make LLMs Generate 100% Valid Molecules?* arXiv:2509.23099.
- Cao N.D., Kipf T. *MolGAN.* arXiv:1805.11973.
- Polykovskiy D. et al. *MOSES.* *Front. Pharmacol.* 2020. https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2020.565644/full. Code: https://github.com/molecularsets/moses.
- Brown N., Fiscato M., Segler M.H.S., Vaucher A.C. *GuacaMol.* *JCIM* 2019. https://pubs.acs.org/doi/10.1021/acs.jcim.8b00839.
- Ruddigkeit L., van Deursen R., Blum L.C., Reymond J.-L. *Enumeration of 166 Billion Organic Small Molecules in GDB-17.* *JCIM* 2012. https://pubs.acs.org/doi/10.1021/ci300415d.
- Krenn M., Häse F., Nigam A., Friederich P., Aspuru-Guzik A. *SELFIES: A 100% robust molecular string representation.* *Mach. Learn.: Sci. Tech.* 1:045024, 2020. arXiv:1905.13741.
- Rong Y. et al. *GROVER.* NeurIPS 2020. arXiv:2007.02835.
- Hu W., Liu B., Gomes J., Zitnik M., Liang P., Pande V., Leskovec J. *Strategies for Pre-training GNNs.* ICLR 2020. arXiv:1905.12265.
- Wellawatte G.P., Seshadri A., White A.D. *MMACE: model-agnostic counterfactual explanations for molecules.* *Chem. Sci.* 13:3697, 2022. https://doi.org/10.1039/D1SC05259D.
- Dalke A., Hert J., Kramer C. *mmpdb.* *JCIM* 2018. https://pubs.acs.org/doi/10.1021/acs.jcim.8b00173.
- Acosta Murillo F. et al. *Benchmarking Molecular Mutation Operators for Evolutionary Drug Design.* *IJMS* 2025. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12692386/.
- Tan H. et al. *NoiseMol.* *J. Mol. Graph. Mod.* 2023. https://www.sciencedirect.com/science/article/abs/pii/S1093326323000529.
- Cooper J.B. et al. *A quantitative definition of hypervalency.* *Chem. Sci.* 2015. https://pubs.rsc.org/en/content/articlehtml/2015/sc/c5sc02076j.
- SCAGE authors. *Self-conformation-aware pre-training framework.* *Nat. Commun.* 2025. https://www.nature.com/articles/s41467-025-59634-0.
- Bilodeau C., Jin W., Jaakkola T., Barzilay R., Jensen K.F. *Generative Models for Molecular Discovery.* *WIREs Comput. Mol. Sci.* 12:e1608, 2022. https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1608.
