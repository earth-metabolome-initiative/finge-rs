//! Sparse (unfolded) representations for the Tanimoto index.
//!
//! [`SparseFingerprint<T>`] is a sorted, deduplicated set of feature ids, the
//! exact unfolded counterpart of a folded [`crate::BitFingerprint`]. It plugs
//! into [`TanimotoIndex`](super::TanimotoIndex) through the same
//! [`TanimotoItem`] trait, so folded and unfolded fingerprints query through one
//! index. `SparseCountFingerprint` and its impls land in Phase 4.

use alloc::vec::Vec;

use super::TanimotoItem;

/// A sparse fingerprint: a sorted, deduplicated set of feature ids.
///
/// This is the exact unfolded feature set (raw Morgan identifiers for ECFP,
/// shingle keys for MAP4), so its [`TanimotoItem::tanimoto`] is the true
/// Jaccard of the two sets with no folding collisions. The pruning magnitude is
/// the set length.
///
/// ```
/// use finge_rs::{SparseFingerprint, TanimotoItem};
///
/// let a = SparseFingerprint::from_ids(vec![3u32, 1, 2, 1]);
/// let b = SparseFingerprint::from_ids(vec![2u32, 3, 4]);
/// // The builder sorts and deduplicates.
/// assert_eq!(a.as_slice(), &[1, 2, 3]);
/// // Intersection {2, 3} = 2 over union {1, 2, 3, 4} = 4.
/// assert_eq!(a.tanimoto(&b), 0.5);
/// assert_eq!(a.tanimoto(&a), 1.0);
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseFingerprint<T> {
    ids: Vec<T>,
}

impl<T: Ord + Copy> SparseFingerprint<T> {
    /// Builds a sparse fingerprint from a feature-id list, sorting and
    /// deduplicating it so the set invariant holds.
    #[must_use]
    pub fn from_ids(mut ids: Vec<T>) -> Self {
        ids.sort_unstable();
        ids.dedup();
        Self { ids }
    }

    /// Returns the sorted, deduplicated feature ids.
    #[inline]
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.ids
    }

    /// Returns the number of distinct features.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Returns whether the fingerprint carries no features.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }
}

impl<T: Ord + Copy> FromIterator<T> for SparseFingerprint<T> {
    #[inline]
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self::from_ids(iter.into_iter().collect())
    }
}

impl<T: Ord + Copy> TanimotoItem for SparseFingerprint<T> {
    #[inline]
    fn magnitude(&self) -> u64 {
        u64::try_from(self.ids.len()).expect("a set length fits in u64")
    }

    fn tanimoto(&self, other: &Self) -> f64 {
        let intersection = intersection_count(&self.ids, &other.ids);
        let union = self.magnitude() + other.magnitude() - intersection;
        if union == 0 {
            return 0.0;
        }
        // Set cardinalities are far below f64's 2^53 exact-integer ceiling, so
        // the widening ratio is lossless.
        debug_assert!(intersection < (1u64 << 53) && union < (1u64 << 53));
        intersection as f64 / union as f64
    }
}

/// Counts the shared ids of two sorted, deduplicated slices with a galloping
/// (exponential then binary) search. Probing the smaller set against the larger
/// costs `O(|small| * log |large|)`, beating a linear merge when the sizes are
/// very unequal, which is common for unfolded fingerprints.
fn intersection_count<T: Ord + Copy>(a: &[T], b: &[T]) -> u64 {
    let (small, mut large) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let mut count = 0u64;
    for &x in small {
        if large.is_empty() {
            break;
        }
        // Exponential search for the first index whose id is >= x, so the
        // window narrows to `[bound / 2, bound)` before the binary step.
        let mut bound = 1usize;
        while bound < large.len() && large[bound] < x {
            bound <<= 1;
        }
        let lo = bound >> 1;
        let hi = bound.min(large.len());
        let pos = lo + large[lo..hi].partition_point(|&v| v < x);
        // Everything before `pos` is strictly less than x and never matches a
        // later (larger) probe, so discard it. Matches are unique by dedup.
        large = &large[pos..];
        if large.first() == Some(&x) {
            count += 1;
            large = &large[1..];
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use alloc::collections::BTreeSet;
    use alloc::vec::Vec;

    use minhash_rs::prelude::{MinHash, MinHasher};
    use proptest::prelude::*;
    use smiles_parser::smiles::Smiles;

    use super::SparseFingerprint;
    use crate::smiles_support::SmilesRdkitScratch;
    use crate::tanimoto::card_bound;
    use crate::test_fixtures::rdkit_ecfp_fixture;
    use crate::{EcfpFingerprint, LshIndex, TanimotoIndex, TanimotoItem};

    fn sf(ids: &[u32]) -> SparseFingerprint<u32> {
        SparseFingerprint::from_ids(ids.to_vec())
    }

    /// Independent Jaccard reference over `BTreeSet`s, unrelated to the
    /// galloping-merge kernel under test.
    fn set_tanimoto(a: &BTreeSet<u32>, b: &BTreeSet<u32>) -> f64 {
        if a.is_empty() && b.is_empty() {
            return 0.0;
        }
        let intersection = a.intersection(b).count();
        let union = a.len() + b.len() - intersection;
        intersection as f64 / union as f64
    }

    #[test]
    fn tanimoto_matches_hand_values() {
        let a = sf(&[0, 1, 2]);
        let b = sf(&[1, 2, 3]);
        // Intersection {1, 2} = 2 over union {0, 1, 2, 3} = 4.
        assert_eq!(a.tanimoto(&b), 0.5);
        assert_eq!(a.tanimoto(&a), 1.0);
        // Disjoint sets share nothing.
        assert_eq!(sf(&[0, 1]).tanimoto(&sf(&[2, 3])), 0.0);
    }

    #[test]
    fn tanimoto_empty_conventions() {
        let empty = SparseFingerprint::<u32>::from_ids(Vec::new());
        let a = sf(&[0, 1]);
        // Both empty is 0.0 by contract, and an empty side shares nothing.
        assert_eq!(empty.tanimoto(&empty), 0.0);
        assert_eq!(a.tanimoto(&empty), 0.0);
        assert_eq!(empty.tanimoto(&a), 0.0);
    }

    fn arb_set() -> impl Strategy<Value = Vec<u32>> {
        prop::collection::vec(0u32..50, 0..40)
    }

    proptest! {
        /// The galloping-merge Tanimoto matches the independent set reference
        /// exactly across random id lists (with duplicates), including empties
        /// and wildly unequal sizes.
        #[test]
        fn tanimoto_matches_set_reference(a in arb_set(), b in arb_set()) {
            let fa = SparseFingerprint::from_ids(a.clone());
            let fb = SparseFingerprint::from_ids(b.clone());
            let sa: BTreeSet<u32> = a.into_iter().collect();
            let sb: BTreeSet<u32> = b.into_iter().collect();
            prop_assert_eq!(fa.tanimoto(&fb), set_tanimoto(&sa, &sb));
        }

        /// The mandatory safety bound: `card_bound` never underestimates the
        /// exact Tanimoto for the sparse representation.
        #[test]
        fn card_bound_upper_bounds_tanimoto(a in arb_set(), b in arb_set()) {
            let fa = SparseFingerprint::from_ids(a);
            let fb = SparseFingerprint::from_ids(b);
            prop_assert!(card_bound(fa.magnitude(), fb.magnitude()) >= fa.tanimoto(&fb));
        }
    }

    /// Exact unfolded-ECFP Tanimoto tracks the MinHash `estimate_jaccard_index`
    /// within MinHash variance at 512 permutations, confirming the unfolded set
    /// construction represents the same features the sketch samples.
    #[test]
    fn unfolded_ecfp_tanimoto_tracks_minhash() {
        let fixture = rdkit_ecfp_fixture();
        let fingerprint = EcfpFingerprint::new(2, 2048);

        let mut sparses: Vec<SparseFingerprint<u32>> = Vec::new();
        let mut sketches: Vec<MinHash<u32, 512>> = Vec::new();
        for smiles in fixture.molecules.iter().take(256) {
            let mol: Smiles = smiles.parse().expect("fixture SMILES should parse");
            let mut scratch = SmilesRdkitScratch::default();
            let graph = scratch.prepare(&mol);
            sparses.push(fingerprint.unfolded(&graph));
            sketches.push(fingerprint.minhash::<_, u32, 512>(&graph));
        }

        let count = sparses.len();
        let mut pairs = 0usize;
        let mut total_error = 0.0_f64;
        for &stride in &[1usize, 7, 53] {
            for i in 0..count {
                let j = i + stride;
                if j >= count {
                    break;
                }
                let estimate = sketches[i].estimate_jaccard_index(&sketches[j]);
                let exact = sparses[i].tanimoto(&sparses[j]);
                total_error += (estimate - exact).abs();
                pairs += 1;
            }
        }

        let mean_error = total_error / pairs as f64;
        assert!(
            mean_error < 0.02,
            "unfolded ECFP Tanimoto strayed from MinHash estimate: {mean_error}"
        );
    }

    /// A real-corpus guaranteed-k check: the Tanimoto index returns exactly `k`
    /// neighbours (self excluded) for every node, while a strict-threshold LSH
    /// index under-fills at least one neighbourhood, the sparsity failure the
    /// Tanimoto index eliminates.
    #[test]
    fn tanimoto_guarantees_k_where_lsh_does_not() {
        const N: usize = 512;
        let fixture = rdkit_ecfp_fixture();
        let fingerprint = EcfpFingerprint::new(2, 2048);
        let k = 10usize;

        let mut index: TanimotoIndex<SparseFingerprint<u32>> = TanimotoIndex::new();
        let mut sketches: Vec<MinHash<u32, N>> = Vec::new();
        for smiles in fixture.molecules.iter().take(128) {
            let mol: Smiles = smiles.parse().expect("fixture SMILES should parse");
            let mut scratch = SmilesRdkitScratch::default();
            let graph = scratch.prepare(&mol);
            index.insert(fingerprint.unfolded(&graph));
            sketches.push(fingerprint.minhash::<_, u32, N>(&graph));
        }
        index.freeze();

        let n = index.len();
        for id in 0..n {
            let id64 = u64::try_from(id).expect("corpus fits u64");
            let hits = index.query_excluding(index.item(id64), id64, k);
            assert_eq!(hits.len(), k.min(n - 1));
        }

        let mut lsh = LshIndex::<u32, N, 8>::new();
        for sketch in &sketches {
            lsh.insert(*sketch);
        }
        let lsh_under_fills = (0..n).any(|id| {
            let id64 = u64::try_from(id).expect("corpus fits u64");
            lsh.query(&sketches[id], k + 1)
                .into_iter()
                .filter(|(hit, _)| *hit != id64)
                .count()
                < k
        });
        assert!(
            lsh_under_fills,
            "expected the strict LSH index to leave a neighbourhood under-filled"
        );
    }
}
