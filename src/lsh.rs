//! Banded MinHash locality-sensitive hashing (LSH) index.
//!
//! The index is fingerprint-agnostic: it stores and bands MinHash *signatures*,
//! so any fingerprint that can be transformed into a [`MinHash`] plugs in. The
//! transformation is the only fingerprint-specific piece, captured by the
//! [`Sketcher`] trait (identity for an object that is already a signature,
//! shingles-to-MinHash for MAP4, active-bins-to-MinHash for a folded bit
//! fingerprint, and so on).
//!
//! Banding splits the `N` signature registers into `bands` contiguous bands of
//! `N / bands` rows. Two objects are candidate neighbours when they collide in
//! at least one band. The probability that two objects with Jaccard similarity
//! `s` collide is `1 - (1 - s^rows)^bands`, an S-curve whose threshold is tuned
//! by the `bands` and `rows` split. Candidates are then refined by the MinHash
//! Jaccard estimate.

use alloc::{
    collections::{BTreeMap, BTreeSet},
    vec::Vec,
};

use minhash_rs::prelude::{Maximal, MinHash, MinHasher, Primitive};

/// Transformation from an item to a MinHash signature (possibly identity).
pub trait Sketcher<Item, Word, const N: usize>
where
    u64: Primitive<Word>,
{
    /// Returns the MinHash signature for `item`.
    fn sketch(&self, item: &Item) -> MinHash<Word, N>;
}

/// A banded MinHash LSH index over `N`-permutation signatures.
///
/// ```
/// use finge_rs::{LshIndex, MinHash};
///
/// // Two item sets sketched into 128-permutation signatures, and a copy of the
/// // first. The index uses 32 bands of 4 rows each (128 / 32).
/// let first: MinHash<u32, 128> = [1u64, 2, 3, 4, 5].iter().copied().collect();
/// let second: MinHash<u32, 128> = [4u64, 5, 6, 7, 8].iter().copied().collect();
/// let first_again: MinHash<u32, 128> = [1u64, 2, 3, 4, 5].iter().copied().collect();
///
/// let mut index = LshIndex::<u32, 128, 32>::new();
/// index.insert(first);
/// index.insert(second);
/// let query_id = index.insert(first_again);
///
/// // Querying with the copy retrieves the original at estimated similarity 1.0.
/// let hits = index.query(index.signature(query_id), 2);
/// assert_eq!(hits[0].1, 1.0);
/// ```
#[derive(Debug, Clone)]
pub struct LshIndex<Word, const N: usize, const BANDS: usize>
where
    u64: Primitive<Word>,
{
    buckets: Vec<BTreeMap<u64, Vec<u32>>>,
    signatures: Vec<MinHash<Word, N>>,
}

impl<
    Word: Copy + Ord + Maximal + Primitive<u64> + core::hash::Hash,
    const N: usize,
    const BANDS: usize,
> LshIndex<Word, N, BANDS>
where
    u64: Primitive<Word>,
{
    /// Creates an empty index that splits each signature into `BANDS` bands.
    ///
    /// `BANDS` must divide `N`; the divisibility is checked at compile time by
    /// [`MinHash::band_hashes`], the first time a signature is banded.
    #[must_use]
    pub fn new() -> Self {
        Self {
            buckets: (0..BANDS).map(|_| BTreeMap::new()).collect(),
            signatures: Vec::new(),
        }
    }

    /// Returns the number of bands.
    #[inline]
    #[must_use]
    pub const fn bands(&self) -> usize {
        BANDS
    }

    /// Returns the number of indexed signatures.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.signatures.len()
    }

    /// Returns whether the index is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.signatures.is_empty()
    }

    /// Inserts a signature, returning its assigned id.
    pub fn insert(&mut self, signature: MinHash<Word, N>) -> u32 {
        let id = u32::try_from(self.signatures.len()).expect("index holds at most u32::MAX items");
        for (band, key) in signature.band_hashes::<BANDS>().into_iter().enumerate() {
            self.buckets[band].entry(key).or_default().push(id);
        }
        self.signatures.push(signature);
        id
    }

    /// Returns the candidate neighbour ids for `signature`: the union of items
    /// sharing at least one band bucket.
    #[must_use]
    pub fn candidates(&self, signature: &MinHash<Word, N>) -> Vec<u32> {
        let mut candidates = BTreeSet::new();
        for (band, key) in signature.band_hashes::<BANDS>().into_iter().enumerate() {
            if let Some(ids) = self.buckets[band].get(&key) {
                candidates.extend(ids.iter().copied());
            }
        }
        candidates.into_iter().collect()
    }

    /// Returns up to `k` nearest neighbours by MinHash Jaccard estimate among
    /// the banding candidates, best first.
    #[must_use]
    pub fn query(&self, signature: &MinHash<Word, N>, k: usize) -> Vec<(u32, f64)> {
        let mut scored: Vec<(u32, f64)> = self
            .candidates(signature)
            .into_iter()
            .map(|id| {
                (
                    id,
                    self.signatures[id as usize].estimate_jaccard_index(signature),
                )
            })
            .collect();
        scored.sort_unstable_by(|left, right| {
            right
                .1
                .partial_cmp(&left.1)
                .unwrap_or(core::cmp::Ordering::Equal)
                .then(left.0.cmp(&right.0))
        });
        scored.truncate(k);
        scored
    }

    /// Returns the stored signature for an id.
    #[inline]
    #[must_use]
    pub fn signature(&self, id: u32) -> &MinHash<Word, N> {
        &self.signatures[id as usize]
    }
}

impl<
    Word: Copy + Ord + Maximal + Primitive<u64> + core::hash::Hash,
    const N: usize,
    const BANDS: usize,
> Default for LshIndex<Word, N, BANDS>
where
    u64: Primitive<Word>,
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use alloc::{string::String, vec::Vec};
    use std::collections::HashSet;
    use std::time::Instant;

    use minhash_rs::prelude::MinHash;
    use smiles_parser::smiles::Smiles;

    use super::{LshIndex, Sketcher};
    use crate::{fingerprints::Map4Fingerprint, test_fixtures::map4_reference_fixture};

    const N: usize = 512;

    fn jaccard(a: &HashSet<String>, b: &HashSet<String>) -> f64 {
        if a.is_empty() && b.is_empty() {
            return 1.0;
        }
        let intersection = a.iter().filter(|item| b.contains(*item)).count();
        let union = a.len() + b.len() - intersection;
        intersection as f64 / union as f64
    }

    fn corpus() -> (Vec<MinHash<u32, N>>, Vec<HashSet<String>>) {
        let fingerprint = Map4Fingerprint::default();
        let mut sketches = Vec::new();
        let mut sets = Vec::new();
        for mol in &map4_reference_fixture().molecules {
            let parsed: Smiles = mol.smiles.parse().expect("fixture SMILES should parse");
            sketches.push(Sketcher::<_, u32, N>::sketch(&fingerprint, &parsed));
            sets.push(fingerprint.shingles(&parsed).into_iter().collect());
        }
        (sketches, sets)
    }

    fn exact_top_k(sets: &[HashSet<String>], query: usize, k: usize) -> Vec<usize> {
        let mut scored: Vec<(usize, f64)> = (0..sets.len())
            .filter(|&i| i != query)
            .map(|i| (i, jaccard(&sets[query], &sets[i])))
            .collect();
        scored.sort_by(|left, right| {
            right
                .1
                .partial_cmp(&left.1)
                .unwrap()
                .then(left.0.cmp(&right.0))
        });
        scored.truncate(k);
        scored.into_iter().map(|(i, _)| i).collect()
    }

    #[test]
    fn lsh_retrieves_identical_signature_first() {
        let (sketches, _) = corpus();
        let mut index = LshIndex::<u32, N, 64>::new();
        for sketch in &sketches {
            index.insert(*sketch);
        }
        // A molecule queried against the index returns itself, at similarity 1.0.
        for &probe in &[0usize, 5, 100, 500] {
            let hits = index.query(&sketches[probe], 1);
            assert_eq!(hits[0].0 as usize, probe);
            assert_eq!(hits[0].1, 1.0);
        }
    }

    /// Builds an index at a fixed `BANDS`, prints its sweep line, and returns
    /// its recall@k. Factored out because `BANDS` is now a const generic, so the
    /// old runtime `for &bands` loop unrolls into one call per band count.
    fn sweep_bands<const BANDS: usize>(
        sketches: &[MinHash<u32, N>],
        queries: &[usize],
        truths: &[HashSet<usize>],
        k: usize,
        count: usize,
    ) -> f64 {
        let build_start = Instant::now();
        let mut index = LshIndex::<u32, N, BANDS>::new();
        for sketch in sketches {
            index.insert(*sketch);
        }
        let build_time = build_start.elapsed();

        let mut total_recall = 0.0_f64;
        let mut total_candidate_recall = 0.0_f64;
        let mut total_candidates = 0usize;
        let query_start = Instant::now();
        for (&query, truth) in queries.iter().zip(truths) {
            let candidates: HashSet<usize> = index
                .candidates(&sketches[query])
                .into_iter()
                .map(|id| id as usize)
                .collect();
            total_candidates += candidates.len();
            total_candidate_recall += truth.intersection(&candidates).count() as f64 / k as f64;

            let retrieved: HashSet<usize> = index
                .query(&sketches[query], k + 1)
                .into_iter()
                .map(|(id, _)| id as usize)
                .filter(|&id| id != query)
                .take(k)
                .collect();
            total_recall += truth.intersection(&retrieved).count() as f64 / k as f64;
        }
        let query_time = query_start.elapsed();

        let recall = total_recall / queries.len() as f64;
        let candidate_recall = total_candidate_recall / queries.len() as f64;
        let mean_candidates = total_candidates as f64 / queries.len() as f64;
        let rows = N / BANDS;
        std::eprintln!(
            "  bands={BANDS:>3} rows={rows:>2}: recall@{k}={recall:.3} cand_recall={candidate_recall:.3} mean_cand={mean_candidates:>6.1} ({:>4.1}% corpus) build={build_time:?} lsh_query={:.3} ms/q",
            100.0 * mean_candidates / count as f64,
            query_time.as_secs_f64() * 1000.0 / queries.len() as f64,
        );
        recall
    }

    #[test]
    fn lsh_recall_and_selectivity_sweep() {
        let (sketches, sets) = corpus();
        let count = sketches.len();
        let k = 10;
        // A fixed, deterministic query sample.
        let queries: Vec<usize> = (0..count).step_by(10).collect();

        // Precompute exact top-k ground truth once (the brute-force baseline).
        let brute_start = Instant::now();
        let truths: Vec<HashSet<usize>> = queries
            .iter()
            .map(|&query| exact_top_k(&sets, query, k).into_iter().collect())
            .collect();
        let brute_time = brute_start.elapsed();

        std::eprintln!(
            "LSH sweep: N={N} molecules={count} queries={} k={k} | brute-force all-pairs top-{k}: {brute_time:?} ({:.3} ms/query)",
            queries.len(),
            brute_time.as_secs_f64() * 1000.0 / queries.len() as f64,
        );
        let best_recall = sweep_bands::<64>(&sketches, &queries, &truths, k, count)
            .max(sweep_bands::<128>(&sketches, &queries, &truths, k, count))
            .max(sweep_bands::<256>(&sketches, &queries, &truths, k, count))
            .max(sweep_bands::<512>(&sketches, &queries, &truths, k, count));

        assert!(best_recall > 0.5, "no band config reached usable recall");
    }
}

#[cfg(test)]
mod proptests {
    use alloc::vec::Vec;
    use std::collections::HashSet;

    use minhash_rs::prelude::MinHash;
    use proptest::prelude::*;

    use super::LshIndex;

    const N: usize = 64;

    fn sketch(items: &[u64]) -> MinHash<u32, N> {
        items.iter().copied().collect()
    }

    /// Item sets drawn from a small universe so the resulting sketches overlap
    /// across a spread of similarities.
    fn item_sets(min_len: usize) -> impl Strategy<Value = Vec<Vec<u64>>> {
        prop::collection::vec(prop::collection::vec(0u64..40, 0..25), min_len..25)
    }

    /// Fixed band split for the proptests (divides `N`). The recall sweep test
    /// covers the band count as a knob; here it is held constant because it is
    /// now a const generic, and the invariants under test hold for any valid
    /// split.
    const BANDS: usize = 8;

    proptest! {
        /// The index never loses an inserted item from its own bucket, and every
        /// query result is a distinct candidate, ranked by descending similarity
        /// in `[0, 1]`, at most `k` of them.
        #[test]
        fn query_is_well_formed(sets in item_sets(1), k in 1usize..30) {
            let sigs: Vec<MinHash<u32, N>> = sets.iter().map(|s| sketch(s)).collect();
            let mut index = LshIndex::<u32, N, BANDS>::new();
            for sig in &sigs {
                index.insert(*sig);
            }

            for (id, sig) in sigs.iter().enumerate() {
                let candidates = index.candidates(sig);
                prop_assert!(
                    candidates.contains(&(id as u32)),
                    "an inserted item must be a candidate for its own signature"
                );
                let candidate_set: HashSet<u32> = candidates.iter().copied().collect();

                let hits = index.query(sig, k);
                prop_assert!(hits.len() <= k);
                prop_assert!(hits.windows(2).all(|w| w[0].1 >= w[1].1));
                let mut seen = HashSet::new();
                for (hit_id, similarity) in &hits {
                    prop_assert!(candidate_set.contains(hit_id));
                    prop_assert!(seen.insert(*hit_id), "query ids must be distinct");
                    prop_assert!((0.0..=1.0).contains(similarity));
                }
            }
        }
    }
}
