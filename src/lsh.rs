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

use minhash_rs::prelude::MinHash;

/// Transformation from an item to a MinHash signature (possibly identity).
pub trait Sketcher<Item, Word, const N: usize> {
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
/// let first: MinHash<u32, 128> = [1u64, 2, 3, 4, 5].iter().collect();
/// let second: MinHash<u32, 128> = [4u64, 5, 6, 7, 8].iter().collect();
/// let first_again: MinHash<u32, 128> = [1u64, 2, 3, 4, 5].iter().collect();
///
/// let mut index = LshIndex::<u32, 128>::new(32);
/// index.insert(first);
/// index.insert(second);
/// let query_id = index.insert(first_again);
///
/// // Querying with the copy retrieves the original at estimated similarity 1.0.
/// let hits = index.query(index.signature(query_id), 2);
/// assert_eq!(hits[0].1, 1.0);
/// ```
#[derive(Debug, Clone)]
pub struct LshIndex<Word, const N: usize> {
    rows_per_band: usize,
    buckets: Vec<BTreeMap<u64, Vec<u32>>>,
    signatures: Vec<MinHash<Word, N>>,
}

impl<Word: Copy + Eq + Into<u64>, const N: usize> LshIndex<Word, N> {
    /// Creates an index that splits each signature into `bands` bands.
    ///
    /// # Panics
    ///
    /// Panics if `bands` is zero or does not divide `N`.
    #[must_use]
    pub fn new(bands: usize) -> Self {
        assert!(bands > 0, "bands must be positive");
        assert!(N % bands == 0, "bands must divide the permutation count N");
        Self {
            rows_per_band: N / bands,
            buckets: (0..bands).map(|_| BTreeMap::new()).collect(),
            signatures: Vec::new(),
        }
    }

    /// Returns the number of bands.
    #[inline]
    #[must_use]
    pub fn bands(&self) -> usize {
        self.buckets.len()
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
        let registers = signature.as_ref();
        for (band, buckets) in self.buckets.iter_mut().enumerate() {
            let key = band_key(self.rows_per_band, band, registers);
            buckets.entry(key).or_default().push(id);
        }
        self.signatures.push(signature);
        id
    }

    /// Returns the candidate neighbour ids for `signature`: the union of items
    /// sharing at least one band bucket.
    #[must_use]
    pub fn candidates(&self, signature: &MinHash<Word, N>) -> Vec<u32> {
        let registers = signature.as_ref();
        let mut candidates = BTreeSet::new();
        for (band, buckets) in self.buckets.iter().enumerate() {
            if let Some(ids) = buckets.get(&band_key(self.rows_per_band, band, registers)) {
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

#[cfg(feature = "tsne")]
impl<Word: Copy + Eq + Into<u64>, const N: usize> LshIndex<Word, N> {
    /// Builds a fixed-`k` nearest-neighbour graph over every indexed item, in the
    /// layout `bhtsne::tSNE::barnes_hut_with_neighbors` consumes: `row[i]` holds
    /// item `i`'s `k` neighbours by ascending distance (`1 - Jaccard` estimate),
    /// self excluded, every row exactly length `k`.
    ///
    /// LSH candidate sets vary in size, so any item with fewer than `k` banding
    /// candidates falls back to an exact scan over all signatures, guaranteeing
    /// exactly `k` neighbours per row (bhtsne requires equal-length rows).
    ///
    /// # Panics
    ///
    /// Panics if `k == 0` or `k >= self.len()`.
    ///
    /// ```
    /// # #[cfg(feature = "tsne")] {
    /// use finge_rs::{LshIndex, MinHash};
    ///
    /// let mut index = LshIndex::<u32, 128>::new(32);
    /// for items in [vec![1u64, 2, 3], vec![2, 3, 4], vec![7, 8, 9], vec![1, 2, 3, 4]] {
    ///     index.insert(items.iter().collect::<MinHash<u32, 128>>());
    /// }
    ///
    /// // One length-`k` row per item, ready for `bhtsne::tSNE::barnes_hut_with_neighbors`.
    /// let graph = index.knn_graph::<f32>(2);
    /// assert_eq!(graph.len(), 4);
    /// assert!(graph.iter().all(|row| row.len() == 2));
    /// # }
    /// ```
    #[must_use]
    pub fn knn_graph<T>(&self, k: usize) -> Vec<Vec<bhtsne::Neighbor<T>>>
    where
        T: num_traits::NumCast + Copy,
    {
        assert!(k > 0, "k must be positive");
        assert!(
            k < self.signatures.len(),
            "k must be smaller than the number of indexed items"
        );
        (0..self.signatures.len())
            .map(|id| self.knn_row(id, k))
            .collect()
    }

    fn knn_row<T>(&self, id: usize, k: usize) -> Vec<bhtsne::Neighbor<T>>
    where
        T: num_traits::NumCast + Copy,
    {
        let signature = &self.signatures[id];
        let mut scored: Vec<(usize, f64)> = self
            .candidates(signature)
            .into_iter()
            .map(|other| other as usize)
            .filter(|&other| other != id)
            .map(|other| {
                (
                    other,
                    self.signatures[other].estimate_jaccard_index(signature),
                )
            })
            .collect();
        // Banding can under-fill the row; fall back to an exact scan when so.
        if scored.len() < k {
            scored = (0..self.signatures.len())
                .filter(|&other| other != id)
                .map(|other| {
                    (
                        other,
                        self.signatures[other].estimate_jaccard_index(signature),
                    )
                })
                .collect();
        }
        // Ascending distance is descending similarity.
        scored.sort_unstable_by(|left, right| {
            right
                .1
                .partial_cmp(&left.1)
                .unwrap_or(core::cmp::Ordering::Equal)
                .then(left.0.cmp(&right.0))
        });
        scored.truncate(k);
        scored
            .into_iter()
            .map(|(other, similarity)| bhtsne::Neighbor {
                index: other,
                distance: T::from(1.0 - similarity).expect("Jaccard distance is a finite f64"),
            })
            .collect()
    }
}

/// FNV-1a hash of one band's register rows.
fn band_key<Word: Copy + Into<u64>>(rows_per_band: usize, band: usize, registers: &[Word]) -> u64 {
    let start = band * rows_per_band;
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for &register in &registers[start..start + rows_per_band] {
        let value: u64 = register.into();
        for byte in value.to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
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
        let mut index = LshIndex::<u32, N>::new(64);
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
        let mut best_recall = 0.0_f64;
        for &bands in &[64usize, 128, 256, 512] {
            let build_start = Instant::now();
            let mut index = LshIndex::<u32, N>::new(bands);
            for sketch in &sketches {
                index.insert(*sketch);
            }
            let build_time = build_start.elapsed();

            let mut total_recall = 0.0_f64;
            let mut total_candidate_recall = 0.0_f64;
            let mut total_candidates = 0usize;
            let query_start = Instant::now();
            for (&query, truth) in queries.iter().zip(&truths) {
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
            let rows = N / bands;
            std::eprintln!(
                "  bands={bands:>3} rows={rows:>2}: recall@{k}={recall:.3} cand_recall={candidate_recall:.3} mean_cand={mean_candidates:>6.1} ({:>4.1}% corpus) build={build_time:?} lsh_query={:.3} ms/q",
                100.0 * mean_candidates / count as f64,
                query_time.as_secs_f64() * 1000.0 / queries.len() as f64,
            );
            best_recall = best_recall.max(recall);
        }

        assert!(best_recall > 0.5, "no band config reached usable recall");
    }

    #[cfg(feature = "tsne")]
    #[test]
    fn lsh_knn_graph_feeds_bhtsne() {
        let (sketches, _) = corpus();
        let subset = &sketches[..200];
        let mut index = LshIndex::<u32, N>::new(256);
        for sketch in subset {
            index.insert(*sketch);
        }

        let k = 15;
        let knn = index.knn_graph::<f32>(k);

        // Shape and content invariants bhtsne relies on.
        assert_eq!(knn.len(), subset.len());
        for (i, row) in knn.iter().enumerate() {
            assert_eq!(row.len(), k, "every row must have length k");
            assert!(
                row.iter()
                    .all(|nb| nb.index != i && nb.index < subset.len())
            );
            assert!(row.windows(2).all(|w| w[0].distance <= w[1].distance));
            assert!(
                row.iter()
                    .all(|nb| nb.distance >= 0.0 && nb.distance.is_finite())
            );
        }

        // Feed the graph straight into bhtsne (no metric evaluations) and run a
        // short optimization; the embedding must be finite and n * 2 long.
        let data: Vec<u32> = (0..subset.len() as u32).collect();
        let mut tsne = bhtsne::tSNE::<f32, u32, 2>::new(&data);
        tsne.perplexity(k as f32 / 3.0).epochs(50);
        tsne.barnes_hut_with_neighbors(0.5_f32, &knn);
        let embedding = tsne.embedding();
        assert_eq!(embedding.len(), subset.len() * 2);
        assert!(embedding.iter().all(|value| value.is_finite()));
    }
}

#[cfg(test)]
mod proptests {
    use alloc::{vec, vec::Vec};
    use std::collections::HashSet;

    use minhash_rs::prelude::MinHash;
    use proptest::prelude::*;

    use super::LshIndex;

    const N: usize = 64;

    fn sketch(items: &[u64]) -> MinHash<u32, N> {
        items.iter().collect()
    }

    /// Item sets drawn from a small universe so the resulting sketches overlap
    /// across a spread of similarities.
    fn item_sets(min_len: usize) -> impl Strategy<Value = Vec<Vec<u64>>> {
        prop::collection::vec(prop::collection::vec(0u64..40, 0..25), min_len..25)
    }

    /// Band counts that divide `N`.
    fn bands() -> impl Strategy<Value = usize> {
        prop::sample::select(vec![1usize, 2, 4, 8, 16, 32, 64])
    }

    proptest! {
        /// The index never loses an inserted item from its own bucket, and every
        /// query result is a distinct candidate, ranked by descending similarity
        /// in `[0, 1]`, at most `k` of them.
        #[test]
        fn query_is_well_formed(sets in item_sets(1), bands in bands(), k in 1usize..30) {
            let sigs: Vec<MinHash<u32, N>> = sets.iter().map(|s| sketch(s)).collect();
            let mut index = LshIndex::<u32, N>::new(bands);
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

    #[cfg(feature = "tsne")]
    proptest! {
        /// The kNN graph fed to bhtsne must have one row per item, each of exactly
        /// `k` distinct in-range neighbours, self excluded, by ascending finite
        /// distance in `[0, 1]`.
        #[test]
        fn knn_graph_rows_satisfy_the_bhtsne_contract(sets in item_sets(5), bands in bands()) {
            let sigs: Vec<MinHash<u32, N>> = sets.iter().map(|s| sketch(s)).collect();
            let n = sigs.len();
            let mut index = LshIndex::<u32, N>::new(bands);
            for sig in &sigs {
                index.insert(*sig);
            }
            let k = 1 + (n - 1) / 2; // always in [1, n)

            let graph = index.knn_graph::<f32>(k);
            prop_assert_eq!(graph.len(), n);
            for (i, row) in graph.iter().enumerate() {
                prop_assert_eq!(row.len(), k);
                let mut seen = HashSet::new();
                let mut previous = f32::NEG_INFINITY;
                for neighbour in row {
                    prop_assert!(neighbour.index != i && neighbour.index < n);
                    prop_assert!(seen.insert(neighbour.index), "neighbour ids must be distinct");
                    prop_assert!(neighbour.distance >= previous, "distances must ascend");
                    previous = neighbour.distance;
                    prop_assert!(neighbour.distance.is_finite());
                    prop_assert!((0.0..=1.0).contains(&neighbour.distance));
                }
            }
        }
    }
}
