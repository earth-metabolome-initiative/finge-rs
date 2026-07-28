//! Exact Tanimoto k-nearest-neighbour index.
//!
//! One generic [`TanimotoIndex`] drives insert and query over any
//! representation that implements [`TanimotoItem`]. Tanimoto distance `1 - T`
//! is a true metric, so exact k-NN with Baldi cardinality pruning returns
//! exactly `k` neighbours per node and stays exact, which the banded MinHash
//! [`crate::LshIndex`] cannot guarantee in sparse regions of chemical space.

use alloc::collections::BinaryHeap;
use alloc::vec::Vec;
use core::cmp::{Ordering, Reverse};

mod dense;
mod graph;
mod inverted;
mod sparse;

pub use sparse::SparseFingerprint;

/// Converts a push index to a row id. Lossless on the 64-bit hosts this crate
/// targets, enforced by the const assert.
#[inline]
fn id_of(i: usize) -> u64 {
    const { assert!(usize::BITS <= u64::BITS, "usize must fit a u64 row id") };
    // Lossless per the const assert above.
    i as u64
}

/// Converts a row id to a slice index. Lossless on the 64-bit hosts this crate
/// targets, enforced by the const assert.
#[inline]
fn index_of(id: u64) -> usize {
    const { assert!(usize::BITS >= u64::BITS, "row ids need a 64-bit usize") };
    // Lossless per the const assert above.
    id as usize
}

/// An item the [`TanimotoIndex`] can store and compare.
///
/// The representation (dense bits, dense counts, sparse ids, sparse counts)
/// plugs into one index through this trait, mirroring how [`crate::LshIndex`]
/// stays fingerprint-agnostic through [`crate::Sketcher`].
pub trait TanimotoItem {
    /// Pruning magnitude: popcount for bits, L1 sum for counts, set length for
    /// sparse ids. Orders the index and bounds similarity.
    fn magnitude(&self) -> u64;

    /// Exact Tanimoto similarity in `[0, 1]`. Nonempty self-similarity is 1.0,
    /// both-empty is 0.0.
    fn tanimoto(&self, other: &Self) -> f64;
}

/// Exact Tanimoto k-nearest-neighbour index over items of type `F`.
///
/// Items are stored by value alongside their pruning magnitudes. Query results
/// are ordered by descending similarity with ties broken by ascending id,
/// identical to [`crate::LshIndex`]. The index is build-once and query-many:
/// `insert` appends in O(1), and later phases add a magnitude-sorted order used
/// to prune the scan without changing results.
///
/// ```
/// use finge_rs::{BitFingerprint, TanimotoIndex};
///
/// let mut left = BitFingerprint::zeros(8);
/// left.set(0);
/// left.set(1);
/// left.set(2);
/// let mut right = BitFingerprint::zeros(8);
/// right.set(1);
/// right.set(2);
/// right.set(3);
/// let left_again = left.clone();
///
/// let mut index = TanimotoIndex::new();
/// index.insert(left);
/// index.insert(right);
/// let query = index.insert(left_again);
///
/// // The identical fingerprint comes back first at exact similarity 1.0, and
/// // the id tie breaks toward the lower stored id.
/// let hits = index.query(index.item(query), 3);
/// assert_eq!(hits[0], (0, 1.0));
/// assert_eq!(hits[2], (1, 0.5));
/// ```
#[derive(Debug, Clone)]
pub struct TanimotoIndex<F> {
    items: Vec<F>,
    magnitudes: Vec<u64>,
    /// Ids sorted by ascending magnitude, populated by [`Self::freeze`]. Empty
    /// or shorter than `items` (stale after an insert) means a query falls back
    /// to the exact full scan.
    order: Vec<u64>,
}

impl<F> TanimotoIndex<F> {
    /// Creates an empty index.
    #[must_use]
    pub fn new() -> Self {
        Self {
            items: Vec::new(),
            magnitudes: Vec::new(),
            order: Vec::new(),
        }
    }

    /// Returns the number of indexed items.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Returns whether the index is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Returns the stored item for an id.
    #[inline]
    #[must_use]
    pub fn item(&self, id: u64) -> &F {
        &self.items[index_of(id)]
    }
}

impl<F> Default for TanimotoIndex<F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<F: TanimotoItem> TanimotoIndex<F> {
    /// Inserts an item, returning its stable id (its push index). Any prior
    /// freeze is invalidated, so the index falls back to the full scan until it
    /// is frozen again.
    pub fn insert(&mut self, item: F) -> u64 {
        let id = id_of(self.items.len());
        self.magnitudes.push(item.magnitude());
        self.items.push(item);
        self.order.clear();
        id
    }

    /// Precomputes the magnitude-sorted order that the pruned query sweep uses.
    /// Idempotent. An unfrozen or stale index transparently falls back to the
    /// exact full scan, so freezing only changes speed, never results.
    pub fn freeze(&mut self) {
        let mut order: Vec<u64> = (0..self.items.len()).map(id_of).collect();
        order.sort_unstable_by_key(|&id| self.magnitudes[index_of(id)]);
        self.order = order;
    }

    /// Whether a current magnitude-sorted order is available for pruning.
    fn is_frozen(&self) -> bool {
        !self.items.is_empty() && self.order.len() == self.items.len()
    }

    /// Returns up to `min(k, len)` nearest neighbours of `q` by exact Tanimoto
    /// similarity, best first, with ties broken by ascending id.
    #[must_use]
    pub fn query(&self, q: &F, k: usize) -> Vec<(u64, f64)> {
        if self.is_frozen() {
            self.pruned(q, None, k)
        } else {
            self.ranked(q, None, k)
        }
    }

    /// Returns up to `min(k, len - 1)` nearest neighbours of `q`, excluding the
    /// item `skip`, best first, with ties broken by ascending id.
    #[must_use]
    pub fn query_excluding(&self, q: &F, skip: u64, k: usize) -> Vec<(u64, f64)> {
        if self.is_frozen() {
            self.pruned(q, Some(skip), k)
        } else {
            self.ranked(q, Some(skip), k)
        }
    }

    /// Full scan scoring every item against `q`, optionally skipping one id,
    /// ordered by the contract ordering and truncated to `k`. Phase 2 adds the
    /// pruned path that this scan stays the oracle for.
    fn ranked(&self, q: &F, skip: Option<u64>, k: usize) -> Vec<(u64, f64)> {
        // The magnitude column tracks the item column one for one. Reading it
        // here upholds the invariant that the Phase 2 pruning sweep relies on.
        debug_assert_eq!(
            self.magnitudes.len(),
            self.items.len(),
            "items and magnitudes stay in lockstep"
        );
        let mut scored: Vec<(u64, f64)> = self
            .items
            .iter()
            .enumerate()
            .map(|(i, item)| (id_of(i), item))
            .filter(|(id, _)| skip != Some(*id))
            .map(|(id, item)| (id, q.tanimoto(item)))
            .collect();
        scored.sort_unstable_by(|left, right| {
            right
                .1
                .partial_cmp(&left.1)
                .unwrap_or(Ordering::Equal)
                .then(left.0.cmp(&right.0))
        });
        scored.truncate(k);
        scored
    }

    /// Baldi-pruned outward sweep over the magnitude-sorted `order`. Starting at
    /// the query magnitude it walks toward smaller and larger magnitudes,
    /// evaluating the exact `tanimoto` and stopping a direction once
    /// `card_bound` falls below the running k-th-best similarity. The visited
    /// set is a superset of the true top-k, so the shared sort and truncate
    /// yield results identical to the full scan.
    fn pruned(&self, q: &F, skip: Option<u64>, k: usize) -> Vec<(u64, f64)> {
        if k == 0 {
            return Vec::new();
        }
        let mq = q.magnitude();
        let n = self.order.len();
        // First position whose magnitude is at least `mq`. `right` sweeps up
        // from there, `left` sweeps down from just before it.
        let start = self
            .order
            .partition_point(|&id| self.magnitudes[index_of(id)] < mq);
        let mut right = start;
        let mut left = start;
        let mut candidates: Vec<(u64, f64)> = Vec::new();
        // Min-heap of the best k similarities seen. Its smallest element is the
        // running k-th-best score that gates each direction.
        let mut best: BinaryHeap<Reverse<Sim>> = BinaryHeap::new();
        loop {
            let tau = if best.len() == k {
                best.peek().expect("a full heap has a minimum").0.0
            } else {
                0.0
            };
            let up =
                (right < n).then(|| card_bound(mq, self.magnitudes[index_of(self.order[right])]));
            let down =
                (left > 0).then(|| card_bound(mq, self.magnitudes[index_of(self.order[left - 1])]));
            let up_live = up.is_some_and(|bound| bound >= tau);
            let down_live = down.is_some_and(|bound| bound >= tau);
            let go_up = match (up_live, down_live) {
                (true, false) => true,
                (false, true) => false,
                (true, true) => up.unwrap_or(0.0) >= down.unwrap_or(0.0),
                (false, false) => break,
            };
            let pos = if go_up {
                let p = right;
                right += 1;
                p
            } else {
                left -= 1;
                left
            };
            let id = self.order[pos];
            if skip == Some(id) {
                continue;
            }
            let sim = q.tanimoto(&self.items[index_of(id)]);
            candidates.push((id, sim));
            if best.len() < k {
                best.push(Reverse(Sim(sim)));
            } else if let Some(&Reverse(worst)) = best.peek() {
                if Sim(sim) > worst {
                    best.pop();
                    best.push(Reverse(Sim(sim)));
                }
            }
        }
        candidates.sort_unstable_by(|lhs, rhs| {
            rhs.1
                .partial_cmp(&lhs.1)
                .unwrap_or(Ordering::Equal)
                .then(lhs.0.cmp(&rhs.0))
        });
        candidates.truncate(k);
        candidates
    }
}

/// Baldi cardinality upper bound on the Tanimoto similarity of two items given
/// only their magnitudes. Returns 0.0 when either magnitude is zero, otherwise
/// `min(a, b) / max(a, b)`. This is a monotone upper bound for both the binary
/// and the L1-count Tanimoto forms, so one bound function prunes every
/// representation.
#[inline]
fn card_bound(a: u64, b: u64) -> f64 {
    if a == 0 || b == 0 {
        return 0.0;
    }
    let lo = core::cmp::min(a, b);
    let hi = core::cmp::max(a, b);
    // Fingerprint magnitudes sit far below f64's 2^53 exact-integer ceiling, so
    // widening each cardinality to compute the ratio is lossless.
    debug_assert!(
        hi < (1u64 << 53),
        "magnitude exceeds the f64 exact-integer range"
    );
    lo as f64 / hi as f64
}

/// Total order over similarities so a [`BinaryHeap`] can track the running
/// k-th-best score. Similarities are finite values in `[0, 1]`, so `total_cmp`
/// is a full order over them.
#[derive(Clone, Copy, PartialEq)]
struct Sim(f64);

impl Eq for Sim {}

impl PartialOrd for Sim {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Sim {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

#[cfg(test)]
mod tests {
    use alloc::collections::BTreeSet;
    use alloc::vec::Vec;

    use proptest::prelude::*;

    use super::TanimotoIndex;
    use crate::BitFingerprint;

    fn fp(len: usize, bits: &[usize]) -> BitFingerprint {
        let mut fingerprint = BitFingerprint::zeros(len);
        for &bit in bits {
            fingerprint.set(bit);
        }
        fingerprint
    }

    /// Independent Jaccard reference over the set of active bits, deliberately
    /// unrelated to the word-kernel `tanimoto` under test. This is the permanent
    /// oracle that every later phase is checked against.
    fn naive_tanimoto(a: &BitFingerprint, b: &BitFingerprint) -> f64 {
        let sa: BTreeSet<usize> = a.active_bits().collect();
        let sb: BTreeSet<usize> = b.active_bits().collect();
        if sa.is_empty() && sb.is_empty() {
            return 0.0;
        }
        let inter = sa.intersection(&sb).count();
        let union = sa.len() + sb.len() - inter;
        inter as f64 / union as f64
    }

    /// Brute-force top-k by the naive reference, ordered like the index
    /// (descending similarity, ascending id).
    fn brute_knn(
        items: &[BitFingerprint],
        q: &BitFingerprint,
        skip: Option<u64>,
        k: usize,
    ) -> Vec<(u64, f64)> {
        let mut scored: Vec<(u64, f64)> = items
            .iter()
            .enumerate()
            .map(|(i, item)| (u64::try_from(i).expect("corpus fits u64"), item))
            .filter(|(id, _)| skip != Some(*id))
            .map(|(id, item)| (id, naive_tanimoto(q, item)))
            .collect();
        scored.sort_by(|left, right| {
            right
                .1
                .partial_cmp(&left.1)
                .expect("scores are finite")
                .then(left.0.cmp(&right.0))
        });
        scored.truncate(k);
        scored
    }

    #[test]
    fn query_ranks_by_exact_tanimoto() {
        let corpus = [
            fp(8, &[0, 1, 2]),
            fp(8, &[1, 2, 3]),
            fp(8, &[4, 5]),
            fp(8, &[0, 1, 2]),
        ];
        let mut index = TanimotoIndex::new();
        for item in &corpus {
            index.insert(item.clone());
        }
        // A query equal to items 0 and 3 retrieves both at 1.0, ascending id,
        // then item 1 at 0.5 (two shared over a union of four).
        let hits = index.query(&fp(8, &[0, 1, 2]), 3);
        assert_eq!(hits[0], (0, 1.0));
        assert_eq!(hits[1], (3, 1.0));
        assert_eq!(hits[2], (1, 0.5));
    }

    #[test]
    fn query_truncates_and_handles_empty() {
        let mut index = TanimotoIndex::new();
        // An empty index yields nothing.
        assert!(index.query(&fp(8, &[0]), 5).is_empty());
        index.insert(fp(8, &[0, 1]));
        index.insert(fp(8, &[]));
        // A k larger than the corpus returns every item.
        assert_eq!(index.query(&fp(8, &[0, 1]), 10).len(), 2);
        // An empty query scores 0.0 against all, ordered by ascending id.
        assert_eq!(
            index.query(&fp(8, &[]), 10),
            alloc::vec![(0u64, 0.0), (1u64, 0.0)]
        );
    }

    #[test]
    fn query_excluding_drops_the_skipped_id() {
        let corpus = [fp(8, &[0, 1, 2]), fp(8, &[0, 1, 2]), fp(8, &[3])];
        let mut index = TanimotoIndex::new();
        for item in &corpus {
            index.insert(item.clone());
        }
        let hits = index.query_excluding(&fp(8, &[0, 1, 2]), 0, 5);
        assert!(hits.iter().all(|(id, _)| *id != 0));
        assert_eq!(hits[0], (1, 1.0));
    }

    #[test]
    fn freeze_matches_full_scan() {
        let corpus = [
            fp(16, &[0, 1, 2, 9]),
            fp(16, &[1, 2, 3]),
            fp(16, &[4, 5, 6, 7, 8]),
            fp(16, &[0, 1, 2, 9]),
            fp(16, &[]),
        ];
        let mut full = TanimotoIndex::new();
        let mut pruned = TanimotoIndex::new();
        for item in &corpus {
            full.insert(item.clone());
            pruned.insert(item.clone());
        }
        pruned.freeze();
        // Freezing is a pure speedup: every query and every k, including the
        // k > len edge, matches the exact full scan.
        for q in &corpus {
            for k in 0..=corpus.len() + 1 {
                assert_eq!(full.query(q, k), pruned.query(q, k));
                assert_eq!(
                    full.query_excluding(q, 0, k),
                    pruned.query_excluding(q, 0, k)
                );
            }
        }
    }

    #[test]
    fn insert_after_freeze_falls_back_to_full_scan() {
        let mut index = TanimotoIndex::new();
        index.insert(fp(8, &[0, 1]));
        index.freeze();
        // A later insert makes the frozen order stale, so the query must fall
        // back to the exact full scan rather than trust the stale order.
        index.insert(fp(8, &[0, 1, 2]));
        assert_eq!(
            index.query(&fp(8, &[0, 1, 2]), 2),
            alloc::vec![(1u64, 1.0), (0u64, 2.0 / 3.0)]
        );
    }

    fn arb_corpus() -> impl Strategy<Value = (Vec<BitFingerprint>, BitFingerprint)> {
        (1usize..=100).prop_flat_map(|len| {
            let item = prop::collection::vec(0..len, 0..len).prop_map(move |bits| fp(len, &bits));
            (prop::collection::vec(item.clone(), 1..20), item)
        })
    }

    proptest! {
        /// The full-scan `query` matches the independent brute-force oracle
        /// exactly, ids and scores, across random corpora and `k`.
        #[test]
        fn query_equals_brute_force((corpus, q) in arb_corpus(), k in 0usize..25) {
            let mut index = TanimotoIndex::new();
            for item in &corpus {
                index.insert(item.clone());
            }
            prop_assert_eq!(index.query(&q, k), brute_knn(&corpus, &q, None, k));
        }

        /// `query_excluding` matches the oracle with the skipped id removed,
        /// covering the empty-fingerprint and `k > len` edges.
        #[test]
        fn query_excluding_equals_brute_force(
            (corpus, q) in arb_corpus(),
            k in 0usize..25,
            raw_skip in 0usize..64,
        ) {
            let mut index = TanimotoIndex::new();
            for item in &corpus {
                index.insert(item.clone());
            }
        let skip = u64::try_from(raw_skip % corpus.len()).expect("corpus fits u64");
            prop_assert_eq!(
                index.query_excluding(&q, skip, k),
                brute_knn(&corpus, &q, Some(skip), k)
            );
        }

        /// The pruned outward sweep on a frozen index returns byte-for-byte the
        /// same ids and scores as the brute-force oracle, across random corpora,
        /// `k`, and query magnitudes, including empty fingerprints and `k > len`.
        #[test]
        fn pruned_query_equals_brute_force((corpus, q) in arb_corpus(), k in 0usize..25) {
            let mut index = TanimotoIndex::new();
            for item in &corpus {
                index.insert(item.clone());
            }
            index.freeze();
            prop_assert_eq!(index.query(&q, k), brute_knn(&corpus, &q, None, k));
        }

        /// `query_excluding` on a frozen index matches the oracle with the
        /// skipped id removed.
        #[test]
        fn pruned_query_excluding_equals_brute_force(
            (corpus, q) in arb_corpus(),
            k in 0usize..25,
            raw_skip in 0usize..64,
        ) {
            let mut index = TanimotoIndex::new();
            for item in &corpus {
                index.insert(item.clone());
            }
            index.freeze();
        let skip = u64::try_from(raw_skip % corpus.len()).expect("corpus fits u64");
            prop_assert_eq!(
                index.query_excluding(&q, skip, k),
                brute_knn(&corpus, &q, Some(skip), k)
            );
        }
    }
}
