//! MAP4 (MinHashed Atom Pair, diameter 4) fingerprints.
//!
//! MAP4 represents a molecule as a set of shingles, one per (atom pair, radius).
//! Each shingle is `min(CSa, CSb) | distance | max(CSa, CSb)` where `CSr(i)` is
//! the rooted, non-isomeric circular-substructure SMILES of the radius-`r`
//! environment of atom `i` (the empty label for an empty environment), and
//! `distance` is the topological shortest-path distance between the two atoms
//! (RDKit's `1e8` sentinel for disconnected pairs).
//!
//! The deduplicated shingle set is the dialect-dependent core, validated against
//! an RDKit oracle for environment-partition and Tanimoto parity. Three
//! representations are layered on top:
//!
//! - [`Map4Fingerprint`] folds the shingles into a [`BitFingerprint`] (presence
//!   per bin).
//! - [`CountMap4Fingerprint`] folds them into a [`CountFingerprint`] (shingle
//!   multiplicity per bin, like
//!   [`CountEcfpFingerprint`](super::CountEcfpFingerprint)).
//! - [`Map4Fingerprint::minhash`] builds a MinHash sketch whose
//!   `estimate_jaccard_index` approximates the shingle-set Jaccard similarity.

use alloc::{collections::BTreeSet, format, string::String, vec::Vec};

use geometric_traits::traits::{DenseValuedMatrix, MonoplexGraph, algorithms::PairwiseBFS};
use minhash_rs::prelude::{Maximal, MinHash, MinHasher, Primitive};

use crate::{
    bit_fingerprint::BitFingerprint,
    count_fingerprint::CountFingerprint,
    fingerprint::Fingerprint,
    tanimoto::SparseFingerprint,
    traits::{Map4Graph, MolecularAtom, MolecularBond},
};

/// RDKit `GetDistanceMatrix` sentinel for topologically disconnected atom pairs.
pub const MAP4_DISCONNECTED_DISTANCE: usize = 100_000_000;

/// Default MAP4 radius (diameter 4 means radius 2).
pub const DEFAULT_MAP4_RADIUS: usize = 2;

/// Default folded MAP4 length (the reference MAP4 dimension).
pub const DEFAULT_MAP4_FP_SIZE: usize = 1024;

/// Folded bit MAP4 fingerprint, and the shared shingle and MinHash core.
///
/// ```
/// use finge_rs::{Fingerprint, Map4Fingerprint, MinHasher};
/// use finge_rs::smiles_support::SmilesRdkitScratch;
/// use smiles_parser::smiles::Smiles;
///
/// let molecule: Smiles = "OCC1OC(O)C(O)C(O)C1O".parse().unwrap();
/// let mut scratch = SmilesRdkitScratch::default();
/// let graph = scratch.prepare(&molecule);
///
/// let fingerprint = Map4Fingerprint::default();
/// let shingles = fingerprint.shingles(&graph);      // the raw substructure set
/// let _bits = fingerprint.compute(&graph);          // folded into a `BitFingerprint`
/// let signature = fingerprint.minhash::<_, u32, 1024>(&graph); // a MinHash sketch
///
/// assert!(!shingles.is_empty());
/// assert_eq!(signature.estimate_jaccard_index(&signature), 1.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Map4Fingerprint {
    radius: usize,
    fp_size: usize,
}

impl Map4Fingerprint {
    /// Creates a MAP4 fingerprint with the given environment radius and folded
    /// length.
    #[inline]
    #[must_use]
    pub const fn new(radius: usize, fp_size: usize) -> Self {
        Self { radius, fp_size }
    }

    /// Returns the environment radius.
    #[inline]
    #[must_use]
    pub const fn radius(self) -> usize {
        self.radius
    }

    /// Returns the folded bit-vector length.
    #[inline]
    #[must_use]
    pub const fn fp_size(self) -> usize {
        self.fp_size
    }

    /// Produces the deduplicated, lexicographically sorted MAP4 shingle set for
    /// `graph`.
    ///
    /// Mirrors the reference `map4` algorithm with `is_counted = false`:
    /// environments at radii `1..=radius`, unordered atom pairs `i < j`, the two
    /// environment labels sorted before joining, and the whole set deduplicated.
    #[must_use]
    pub fn shingles<G>(&self, graph: &G) -> Vec<String>
    where
        G: Map4Graph<NodeId = usize>,
        G::NodeSymbol: MolecularAtom,
        G::Bond: MolecularBond<NodeId = usize>,
        <G as MonoplexGraph>::Edges: PairwiseBFS,
    {
        let mut shingles: BTreeSet<String> = BTreeSet::new();
        self.for_each_shingle(graph, |shingle| {
            shingles.insert(shingle);
        });
        shingles.into_iter().collect()
    }

    /// Builds a MinHash sketch of this molecule's shingle set. Compare two
    /// sketches with `estimate_jaccard_index` to approximate their MAP4 Tanimoto
    /// similarity.
    ///
    /// `Word` selects the register width (a memory and precision tradeoff) and
    /// `N` the number of permutations, e.g. `fp.minhash::<_, u32, 1024>(&graph)`.
    #[must_use]
    pub fn minhash<G, Word, const N: usize>(&self, graph: &G) -> MinHash<Word, N>
    where
        G: Map4Graph<NodeId = usize>,
        G::NodeSymbol: MolecularAtom,
        G::Bond: MolecularBond<NodeId = usize>,
        <G as MonoplexGraph>::Edges: PairwiseBFS,
        Word: Ord + Copy + Maximal + Primitive<u64>,
        u64: Primitive<Word>,
    {
        // Gather the shingle keys, then deduplicate before inserting: each
        // `insert_with_siphashes13` costs `O(N)` register updates, so folding the
        // `C(n, 2) * radius` occurrences down to distinct keys first is what keeps
        // sketching fast. Integer sort over `u64` keys replaces the old
        // `BTreeSet<String>` dedup, dropping the per-occurrence heap string.
        let mut keys: Vec<u64> = Vec::new();
        self.for_each_shingle_key(graph, |key| keys.push(key));
        keys.sort_unstable();
        keys.dedup();

        let mut sketch = MinHash::new();
        for key in keys {
            sketch.insert(key);
        }
        sketch
    }

    /// Builds the exact unfolded MAP4 shingle-key set as a
    /// [`SparseFingerprint`], for exact Tanimoto retrieval through
    /// [`TanimotoIndex`](crate::TanimotoIndex). Unlike the folded
    /// [`compute`](crate::Fingerprint::compute) bit vector, this set carries no
    /// fold collisions, so its Tanimoto is the true shingle-set Jaccard (up to
    /// negligible 64-bit key collisions).
    ///
    /// ```
    /// use finge_rs::{Map4Fingerprint, TanimotoItem};
    /// use finge_rs::smiles_support::SmilesRdkitScratch;
    /// use smiles_parser::smiles::Smiles;
    ///
    /// let molecule: Smiles = "OCC1OC(O)C(O)C(O)C1O".parse().unwrap();
    /// let mut scratch = SmilesRdkitScratch::default();
    /// let graph = scratch.prepare(&molecule);
    ///
    /// let features = Map4Fingerprint::default().unfolded(&graph);
    /// assert!(!features.is_empty());
    /// assert_eq!(features.tanimoto(&features), 1.0);
    /// ```
    #[must_use]
    pub fn unfolded<G>(&self, graph: &G) -> SparseFingerprint<u64>
    where
        G: Map4Graph<NodeId = usize>,
        G::NodeSymbol: MolecularAtom,
        G::Bond: MolecularBond<NodeId = usize>,
        <G as MonoplexGraph>::Edges: PairwiseBFS,
    {
        let mut keys: Vec<u64> = Vec::new();
        self.for_each_shingle_key(graph, |key| keys.push(key));
        SparseFingerprint::from_ids(keys)
    }

    /// Visits every shingle occurrence (before deduplication) exactly once,
    /// mapping each environment label through `make_label` and each canonically
    /// ordered `(low, distance, high)` triple through `combine`. The string set
    /// view and the hashed sketch and fold paths share this one enumeration, so
    /// they can never drift into representing different features.
    fn for_each_shingle_by<G, L, T>(
        &self,
        graph: &G,
        mut make_label: impl FnMut(Option<String>) -> L,
        combine: impl Fn(&L, usize, &L) -> T,
        mut visit: impl FnMut(T),
    ) where
        G: Map4Graph<NodeId = usize>,
        G::NodeSymbol: MolecularAtom,
        G::Bond: MolecularBond<NodeId = usize>,
        <G as MonoplexGraph>::Edges: PairwiseBFS,
        L: Ord,
    {
        let atom_count = graph.atom_count();
        if self.radius == 0 || atom_count < 2 {
            return;
        }

        // labels[atom][r - 1] is the radius-r environment label representation.
        let mut labels: Vec<Vec<L>> = Vec::with_capacity(atom_count);
        for atom in 0..atom_count {
            let mut row = Vec::with_capacity(self.radius);
            for radius in 1..=self.radius {
                row.push(make_label(graph.map4_environment_label(atom, radius)));
            }
            labels.push(row);
        }

        // Unweighted all-pairs shortest paths via geometric-traits. Unreachable
        // pairs come back as `None`. MAP4 renders them with RDKit's `1e8`
        // disconnected sentinel.
        let distances = graph.edges().pairwise_bfs();

        for i in 0..atom_count {
            for j in (i + 1)..atom_count {
                let distance = distances
                    .value((i, j))
                    .unwrap_or(MAP4_DISCONNECTED_DISTANCE);
                for (a, b) in labels[i].iter().zip(&labels[j]) {
                    let (low, high) = if a <= b { (a, b) } else { (b, a) };
                    visit(combine(low, distance, high));
                }
            }
        }
    }

    /// The string set view: joins each canonically ordered environment-label
    /// pair and its distance as `low|distance|high`. This is the RDKit-oracle
    /// validated MAP4 feature representation, kept for `shingles` and its tests.
    fn for_each_shingle<G, F>(&self, graph: &G, visit: F)
    where
        G: Map4Graph<NodeId = usize>,
        G::NodeSymbol: MolecularAtom,
        G::Bond: MolecularBond<NodeId = usize>,
        <G as MonoplexGraph>::Edges: PairwiseBFS,
        F: FnMut(String),
    {
        self.for_each_shingle_by(
            graph,
            |label| label.unwrap_or_default(),
            |low: &String, distance, high: &String| format!("{low}|{distance}|{high}"),
            visit,
        );
    }

    /// The allocation-light hot path: each shingle occurrence becomes a single
    /// `u64` key derived from the two environment-label hashes and the distance,
    /// with no per-pair heap allocation. Equal shingles hash to equal keys, so
    /// the key-set Jaccard matches the string-set Jaccard up to negligible 64-bit
    /// collisions. Callers that need the distinct-key set (the MinHash sketch)
    /// deduplicate the emitted keys themselves.
    fn for_each_shingle_key<G, F>(&self, graph: &G, visit: F)
    where
        G: Map4Graph<NodeId = usize>,
        G::NodeSymbol: MolecularAtom,
        G::Bond: MolecularBond<NodeId = usize>,
        <G as MonoplexGraph>::Edges: PairwiseBFS,
        F: FnMut(u64),
    {
        self.for_each_shingle_by(
            graph,
            |label| fnv1a(label.unwrap_or_default().as_bytes()),
            |&low, distance, &high| combine_key(low, distance, high),
            visit,
        );
    }

    /// Folds each shingle occurrence into a `fp_size` bin, calling `visit` with
    /// the bin index.
    fn fold<G, F>(&self, graph: &G, mut visit: F)
    where
        G: Map4Graph<NodeId = usize>,
        G::NodeSymbol: MolecularAtom,
        G::Bond: MolecularBond<NodeId = usize>,
        <G as MonoplexGraph>::Edges: PairwiseBFS,
        F: FnMut(usize),
    {
        if self.fp_size == 0 {
            return;
        }
        let fp_size = self.fp_size as u64;
        self.for_each_shingle_key(graph, |key| {
            visit((key % fp_size) as usize);
        });
    }
}

impl Default for Map4Fingerprint {
    #[inline]
    fn default() -> Self {
        Self::new(DEFAULT_MAP4_RADIUS, DEFAULT_MAP4_FP_SIZE)
    }
}

impl<G> Fingerprint<G> for Map4Fingerprint
where
    G: Map4Graph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
    <G as MonoplexGraph>::Edges: PairwiseBFS,
{
    type Output = BitFingerprint;

    fn compute(&self, graph: &G) -> Self::Output {
        let mut fingerprint = BitFingerprint::zeros(self.fp_size);
        self.fold(graph, |index| fingerprint.set(index));
        fingerprint
    }
}

impl<G, Word, const N: usize> crate::lsh::Sketcher<G, Word, N> for Map4Fingerprint
where
    G: Map4Graph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
    <G as MonoplexGraph>::Edges: PairwiseBFS,
    Word: Ord + Copy + Maximal + Primitive<u64>,
    u64: Primitive<Word>,
{
    #[inline]
    fn sketch(&self, item: &G) -> MinHash<Word, N> {
        self.minhash(item)
    }
}

/// Folded counted MAP4 fingerprint: shingle multiplicity per bin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CountMap4Fingerprint(Map4Fingerprint);

impl CountMap4Fingerprint {
    /// Creates a counted MAP4 fingerprint with the given radius and folded
    /// length.
    #[inline]
    #[must_use]
    pub const fn new(radius: usize, fp_size: usize) -> Self {
        Self(Map4Fingerprint::new(radius, fp_size))
    }

    /// Returns the environment radius.
    #[inline]
    #[must_use]
    pub const fn radius(self) -> usize {
        self.0.radius()
    }

    /// Returns the folded vector length.
    #[inline]
    #[must_use]
    pub const fn fp_size(self) -> usize {
        self.0.fp_size()
    }
}

impl Default for CountMap4Fingerprint {
    #[inline]
    fn default() -> Self {
        Self::new(DEFAULT_MAP4_RADIUS, DEFAULT_MAP4_FP_SIZE)
    }
}

impl<G> Fingerprint<G> for CountMap4Fingerprint
where
    G: Map4Graph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
    <G as MonoplexGraph>::Edges: PairwiseBFS,
{
    type Output = CountFingerprint;

    fn compute(&self, graph: &G) -> Self::Output {
        let mut fingerprint = CountFingerprint::zeros(self.0.fp_size());
        self.0.fold(graph, |index| fingerprint.increment(index));
        fingerprint
    }
}

/// 64-bit FNV-1a over bytes, used to hash environment labels and combine
/// shingle keys.
#[inline]
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Combines a canonically ordered pair of environment-label hashes and their
/// topological distance into one shingle key. Hashing the little-endian bytes
/// keeps the key stable and platform-independent, and gives each distinct
/// `(low, distance, high)` triple its own key up to 64-bit collisions.
#[inline]
fn combine_key(low: u64, distance: usize, high: u64) -> u64 {
    let mut bytes = [0u8; 24];
    bytes[0..8].copy_from_slice(&low.to_le_bytes());
    bytes[8..16].copy_from_slice(&(distance as u64).to_le_bytes());
    bytes[16..24].copy_from_slice(&high.to_le_bytes());
    fnv1a(&bytes)
}

#[cfg(test)]
mod tests {
    use alloc::{
        collections::{BTreeMap, BTreeSet},
        string::{String, ToString},
        vec,
        vec::Vec,
    };
    use std::collections::HashSet;

    use geometric_traits::traits::{DenseValuedMatrix, MonoplexGraph, algorithms::PairwiseBFS};
    use minhash_rs::prelude::{MinHash, MinHasher};
    use smiles_parser::smiles::Smiles;

    use super::{CountMap4Fingerprint, MAP4_DISCONNECTED_DISTANCE, Map4Fingerprint};
    use crate::{
        fingerprint::Fingerprint, smiles_support_impl::SmilesRdkitScratch,
        test_fixtures::map4_reference_fixture, traits::Map4Graph,
    };

    fn sorted_bytes(label: &str) -> Vec<u8> {
        let mut bytes = label.as_bytes().to_vec();
        bytes.sort_unstable();
        bytes
    }

    #[test]
    fn map4_basic_shingles_match_ethanol() {
        // finge-rs env labels for CCO are byte-identical to RDKit here, so the
        // shingle set is too. Atom 1 (central C) has an empty radius-2 shell.
        let mol: Smiles = "CCO".parse().expect("parse");
        let shingles = Map4Fingerprint::default().shingles(&mol);
        assert_eq!(
            shingles,
            vec![
                "C(C)O|1|CC".to_string(),
                "C(C)O|1|OC".to_string(),
                "CCO|2|OCC".to_string(),
                "CC|2|OC".to_string(),
                "|1|CCO".to_string(),
                "|1|OCC".to_string(),
            ]
        );
    }

    /// I3: finge-rs env labels induce a partition that is a clean coarsening of
    /// RDKit's. It must never split an environment RDKit merges (0 consistency
    /// violations), and any merge must be a pure branch-order difference in
    /// RDKit's order-dependent rooted SMILES (the merged RDKit labels are
    /// anagrams of each other).
    #[test]
    fn map4_environment_partition_parity() {
        let fixture = map4_reference_fixture();

        let mut forward: BTreeMap<String, String> = BTreeMap::new();
        let mut inverse: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        let mut consistency_violations = 0usize;
        let mut atom_count_mismatch = 0usize;

        for mol in &fixture.molecules {
            let parsed: Smiles = mol.smiles.parse().expect("fixture SMILES should parse");
            if parsed.nodes().len() != mol.elements.len() {
                atom_count_mismatch += 1;
                continue;
            }
            for (atom, radii) in mol.envs.iter().enumerate() {
                for (radius_index, rdkit_label) in radii.iter().enumerate() {
                    let finge_label = parsed
                        .map4_environment_label(atom, radius_index + 1)
                        .unwrap_or_default();
                    match forward.get(rdkit_label) {
                        None => {
                            forward.insert(rdkit_label.clone(), finge_label.clone());
                        }
                        Some(existing) if existing != &finge_label => {
                            consistency_violations += 1;
                        }
                        _ => {}
                    }
                    inverse
                        .entry(finge_label)
                        .or_default()
                        .insert(rdkit_label.clone());
                }
            }
        }

        assert_eq!(
            atom_count_mismatch, 0,
            "atom indexing diverged from the oracle"
        );
        assert_eq!(
            consistency_violations, 0,
            "finge-rs split an environment RDKit treats as one label",
        );

        // Every merge must be a pure branch-order (anagram) difference.
        let mut merged_finge_labels = 0usize;
        for (finge_label, rdkit_labels) in &inverse {
            if rdkit_labels.len() < 2 {
                continue;
            }
            merged_finge_labels += 1;
            let canonical = sorted_bytes(rdkit_labels.iter().next().unwrap());
            for rdkit_label in rdkit_labels {
                assert_eq!(
                    sorted_bytes(rdkit_label),
                    canonical,
                    "finge label {finge_label:?} merged non-anagram RDKit labels {rdkit_labels:?}",
                );
            }
        }

        // Clean coarsening: at least one such RDKit order-dependence merge
        // exists in this corpus, and finge has no more labels than RDKit.
        assert!(merged_finge_labels > 0);
        assert!(inverse.len() <= forward.len());
    }

    /// I1: finge-rs all-pairs BFS distances (with the 1e8 disconnected sentinel)
    /// match RDKit's distance matrix exactly.
    #[test]
    fn map4_distance_parity() {
        let fixture = map4_reference_fixture();
        for mol in &fixture.molecules {
            let parsed: Smiles = mol.smiles.parse().expect("fixture SMILES should parse");
            if parsed.nodes().len() != mol.elements.len() {
                continue;
            }
            let distances = parsed.edges().pairwise_bfs();
            for (i, row) in mol.distances.iter().enumerate() {
                for (j, &expected) in row.iter().enumerate() {
                    let observed = distances
                        .value((i, j))
                        .unwrap_or(MAP4_DISCONNECTED_DISTANCE);
                    assert_eq!(
                        i64::try_from(observed).unwrap(),
                        expected,
                        "distance mismatch for {} at ({i}, {j})",
                        mol.smiles,
                    );
                }
            }
        }
    }

    fn jaccard<T: Eq + core::hash::Hash>(a: &HashSet<T>, b: &HashSet<T>) -> f64 {
        if a.is_empty() && b.is_empty() {
            return 1.0;
        }
        let intersection = a.iter().filter(|item| b.contains(*item)).count();
        let union = a.len() + b.len() - intersection;
        intersection as f64 / union as f64
    }

    /// I4: finge-rs pairwise Tanimoto tracks RDKit MAP4 Tanimoto. They are not
    /// identical (RDKit's order-dependence over-splits some environments, so
    /// finge-rs reads slightly more similar), but the deviation is small.
    #[test]
    fn map4_shingle_tanimoto_parity() {
        let fixture = map4_reference_fixture();
        let fingerprint = Map4Fingerprint::default();

        let mut finge_sets: Vec<HashSet<String>> = Vec::with_capacity(fixture.molecules.len());
        let mut rdkit_sets: Vec<HashSet<String>> = Vec::with_capacity(fixture.molecules.len());
        for mol in &fixture.molecules {
            let parsed: Smiles = mol.smiles.parse().expect("fixture SMILES should parse");
            finge_sets.push(fingerprint.shingles(&parsed).into_iter().collect());
            rdkit_sets.push(mol.shingles.iter().cloned().collect());
        }

        let count = finge_sets.len();
        let mut pairs = 0usize;
        let mut exact = 0usize;
        let mut max_deviation = 0.0_f64;
        let mut total_deviation = 0.0_f64;
        for &stride in &[1usize, 7, 53, 211] {
            for i in 0..count {
                let j = i + stride;
                if j >= count {
                    break;
                }
                let deviation = (jaccard(&finge_sets[i], &finge_sets[j])
                    - jaccard(&rdkit_sets[i], &rdkit_sets[j]))
                .abs();
                pairs += 1;
                if deviation == 0.0 {
                    exact += 1;
                }
                total_deviation += deviation;
                max_deviation = max_deviation.max(deviation);
            }
        }

        std::eprintln!(
            "MAP4 Tanimoto parity: pairs={pairs} exact={exact} ({:.1}%) mean_dev={:.6} max_dev={:.6}",
            100.0 * exact as f64 / pairs as f64,
            total_deviation / pairs as f64,
            max_deviation,
        );

        assert!(
            max_deviation < 0.05,
            "MAP4 Tanimoto deviation from RDKit too large: {max_deviation}",
        );
    }

    /// The production path (RDKit-normalized graph) yields exactly the same
    /// shingle set as the raw input-order graph the oracle validates, so MAP4 is
    /// representation-independent. This holds only after the upstream fix that
    /// makes canonicalization ignore the kekule order retained on aromatic bonds
    /// (see docs/upstream-issues/smiles-parser-aromatic-bond-order-fragment.md).
    /// Before it, broken-ring environment fragments rendered differently between
    /// the two equivalent representations.
    #[test]
    fn map4_normalized_path_matches_raw() {
        let fixture = map4_reference_fixture();
        let fingerprint = Map4Fingerprint::default();
        let mut scratch = SmilesRdkitScratch::default();

        for mol in &fixture.molecules {
            let parsed: Smiles = mol.smiles.parse().expect("fixture SMILES should parse");
            let raw: BTreeSet<String> = fingerprint.shingles(&parsed).into_iter().collect();
            let prepared = scratch.prepare(&parsed);
            let normalized: BTreeSet<String> =
                fingerprint.shingles(&prepared).into_iter().collect();
            assert_eq!(
                raw, normalized,
                "normalized path diverged for {}",
                mol.smiles
            );
        }
    }

    #[test]
    fn map4_bit_fingerprint_basics() {
        let mol: Smiles = "OCC1OC(O)C(O)C(O)C1O".parse().expect("parse");
        let fingerprint = Map4Fingerprint::new(2, 2048);
        let bits = fingerprint.compute(&mol);
        assert_eq!(bits.len(), 2048);

        let active: Vec<usize> = bits.active_bits().collect();
        assert!(active.windows(2).all(|window| window[0] < window[1]));
        assert!(active.iter().all(|&bit| bit < 2048));
        // Deterministic.
        assert_eq!(
            active,
            fingerprint.compute(&mol).active_bits().collect::<Vec<_>>()
        );
        // At most one bit per distinct shingle.
        assert!(active.len() <= fingerprint.shingles(&mol).len());
    }

    #[test]
    fn map4_count_fingerprint_basics() {
        let mol: Smiles = "OCC1OC(O)C(O)C(O)C1O".parse().expect("parse");
        let bits = Map4Fingerprint::new(2, 2048).compute(&mol);
        let counts = CountMap4Fingerprint::new(2, 2048).compute(&mol);
        assert_eq!(counts.len(), 2048);

        // Presence agrees with the bit fingerprint.
        let bit_set: HashSet<usize> = bits.active_bits().collect();
        let count_set: HashSet<usize> = counts.active_counts().map(|(index, _)| index).collect();
        assert_eq!(bit_set, count_set);

        // Counts are shingle multiplicity, so they sum to the pre-dedup shingle
        // occurrence count: one per (atom pair, radius), i.e. C(n, 2) * radius.
        let atom_count = mol.nodes().len();
        let expected = atom_count * (atom_count - 1) / 2 * 2;
        let total: u32 = counts.active_counts().map(|(_, count)| count).sum();
        assert_eq!(u64::from(total), expected as u64);
    }

    fn bit_set(smiles: &str, fingerprint: Map4Fingerprint) -> HashSet<usize> {
        let mol: Smiles = smiles.parse().expect("fixture SMILES should parse");
        fingerprint.compute(&mol).active_bits().collect()
    }

    /// Tanimoto over the folded bit fingerprint tracks the exact shingle-set
    /// Jaccard (which itself tracks RDKit MAP4 to 99.6%), within folding noise.
    #[test]
    fn map4_folded_tanimoto_tracks_jaccard() {
        let fixture = map4_reference_fixture();
        let fingerprint = Map4Fingerprint::new(2, 4096);

        let mut bits: Vec<HashSet<usize>> = Vec::with_capacity(fixture.molecules.len());
        let mut sets: Vec<HashSet<String>> = Vec::with_capacity(fixture.molecules.len());
        for mol in &fixture.molecules {
            bits.push(bit_set(&mol.smiles, fingerprint));
            let parsed: Smiles = mol.smiles.parse().expect("fixture SMILES should parse");
            sets.push(fingerprint.shingles(&parsed).into_iter().collect());
        }

        let (pairs, mean_error, max_error) = compare_pairs(&bits, &sets);
        std::eprintln!(
            "MAP4 folded vs exact: pairs={pairs} mean_err={mean_error:.6} max_err={max_error:.6}"
        );
        assert!(
            mean_error < 0.01,
            "folded Tanimoto strayed from exact Jaccard"
        );
    }

    /// MinHash `estimate_jaccard_index` tracks the exact shingle-set Jaccard
    /// within MinHash variance at 1024 permutations.
    #[test]
    fn map4_minhash_estimates_jaccard() {
        let fixture = map4_reference_fixture();
        let fingerprint = Map4Fingerprint::default();

        let mut sketches: Vec<MinHash<u32, 1024>> = Vec::with_capacity(fixture.molecules.len());
        let mut sets: Vec<HashSet<String>> = Vec::with_capacity(fixture.molecules.len());
        for mol in &fixture.molecules {
            let parsed: Smiles = mol.smiles.parse().expect("fixture SMILES should parse");
            sketches.push(fingerprint.minhash::<_, u32, 1024>(&parsed));
            sets.push(fingerprint.shingles(&parsed).into_iter().collect());
        }

        let count = sets.len();
        let mut pairs = 0usize;
        let mut total_error = 0.0_f64;
        let mut max_error = 0.0_f64;
        for &stride in &[1usize, 7, 53] {
            for i in 0..count {
                let j = i + stride;
                if j >= count {
                    break;
                }
                let estimate = sketches[i].estimate_jaccard_index(&sketches[j]);
                let exact = jaccard(&sets[i], &sets[j]);
                let error = (estimate - exact).abs();
                pairs += 1;
                total_error += error;
                max_error = max_error.max(error);
            }
        }

        let mean_error = total_error / pairs as f64;
        std::eprintln!(
            "MAP4 minhash vs exact: pairs={pairs} mean_err={mean_error:.6} max_err={max_error:.6}"
        );
        assert!(
            mean_error < 0.02,
            "MinHash estimate strayed from exact Jaccard"
        );
    }

    #[test]
    fn map4_minhash_self_identity_and_determinism() {
        let mol: Smiles = "OCC1OC(O)C(O)C(O)C1O".parse().expect("parse");
        let fingerprint = Map4Fingerprint::default();
        let first = fingerprint.minhash::<_, u32, 512>(&mol);
        let second = fingerprint.minhash::<_, u32, 512>(&mol);
        assert_eq!(first, second);
        assert_eq!(first.estimate_jaccard_index(&second), 1.0);
    }

    /// Compares folded-bit Tanimoto against exact shingle Jaccard over a fixed
    /// deterministic sample of molecule pairs.
    fn compare_pairs(bits: &[HashSet<usize>], sets: &[HashSet<String>]) -> (usize, f64, f64) {
        let count = sets.len();
        let mut pairs = 0usize;
        let mut total_error = 0.0_f64;
        let mut max_error = 0.0_f64;
        for &stride in &[1usize, 7, 53] {
            for i in 0..count {
                let j = i + stride;
                if j >= count {
                    break;
                }
                let folded = jaccard(&bits[i], &bits[j]);
                let exact = jaccard(&sets[i], &sets[j]);
                let error = (folded - exact).abs();
                pairs += 1;
                total_error += error;
                max_error = max_error.max(error);
            }
        }
        (pairs, total_error / pairs as f64, max_error)
    }
}
