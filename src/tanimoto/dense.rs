//! Dense representations for the Tanimoto index.
//!
//! [`TanimotoItem`](super::TanimotoItem) impls for [`crate::BitFingerprint`]
//! (Phase 1, popcount kernel) and [`crate::CountFingerprint`] (Phase 4, min/max
//! kernel) land here.

use super::TanimotoItem;
use crate::BitFingerprint;

impl TanimotoItem for BitFingerprint {
    #[inline]
    fn magnitude(&self) -> u64 {
        u64::try_from(self.as_bitslice().count_ones()).expect("a popcount fits in u64")
    }

    fn tanimoto(&self, other: &Self) -> f64 {
        // Intersection popcount over the word-aligned overlap. `domain` masks
        // the partial tail word, so dead bits never leak into the count, and
        // zipping stops at the shorter operand, so its surplus words feed only
        // the union through the per-side magnitudes.
        let intersection: u64 = self
            .as_bitslice()
            .domain()
            .zip(other.as_bitslice().domain())
            .map(|(left, right)| u64::from((left & right).count_ones()))
            .sum();
        let union = self.magnitude() + other.magnitude() - intersection;
        if union == 0 {
            // Both fingerprints are empty. The contract fixes this at 0.0.
            return 0.0;
        }
        // Cardinalities are popcounts, far below f64's 2^53 exact-integer
        // ceiling, so the widening ratio is lossless.
        debug_assert!(
            union < (1u64 << 53),
            "cardinality exceeds the f64 exact-integer range"
        );
        intersection as f64 / union as f64
    }
}
#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::super::{TanimotoItem, card_bound};
    use crate::BitFingerprint;

    fn fp(len: usize, bits: &[usize]) -> BitFingerprint {
        let mut fingerprint = BitFingerprint::zeros(len);
        for &bit in bits {
            fingerprint.set(bit);
        }
        fingerprint
    }

    #[test]
    fn magnitude_counts_set_bits() {
        assert_eq!(fp(64, &[0, 5, 63]).magnitude(), 3);
        assert_eq!(fp(64, &[]).magnitude(), 0);
        assert_eq!(BitFingerprint::zeros(0).magnitude(), 0);
    }

    #[test]
    fn tanimoto_matches_hand_values() {
        // Overlap of 2 over a union of 4.
        assert_eq!(fp(8, &[0, 1, 2]).tanimoto(&fp(8, &[1, 2, 3])), 0.5);
        // Identical nonempty sets.
        assert_eq!(fp(8, &[0, 1, 2]).tanimoto(&fp(8, &[0, 1, 2])), 1.0);
        // Disjoint sets.
        assert_eq!(fp(8, &[0, 1]).tanimoto(&fp(8, &[2, 3])), 0.0);
        // Subset: 2 shared over a union of 4.
        assert_eq!(fp(8, &[0, 1]).tanimoto(&fp(8, &[0, 1, 2, 3])), 0.5);
    }

    #[test]
    fn tanimoto_empty_conventions() {
        // Both empty is 0.0 by contract, not 1.0.
        assert_eq!(fp(8, &[]).tanimoto(&fp(8, &[])), 0.0);
        // One empty is 0.0.
        assert_eq!(fp(8, &[]).tanimoto(&fp(8, &[0])), 0.0);
        assert_eq!(fp(8, &[0]).tanimoto(&fp(8, &[])), 0.0);
    }

    #[test]
    fn tanimoto_handles_unequal_lengths() {
        // The same set stored at different bit lengths compares as identical.
        assert_eq!(fp(8, &[0, 1]).tanimoto(&fp(128, &[0, 1])), 1.0);
        // Bits beyond the shorter length only enlarge the union.
        assert_eq!(fp(8, &[0, 1, 7]).tanimoto(&fp(128, &[0, 1, 7, 64])), 0.75);
    }

    fn arb_fp() -> impl Strategy<Value = BitFingerprint> {
        (1usize..=128).prop_flat_map(|len| {
            prop::collection::vec(0..len, 0..len).prop_map(move |bits| fp(len, &bits))
        })
    }

    proptest! {
        /// The Baldi bound never underestimates the exact Tanimoto, the safety
        /// property that keeps pruning from dropping true neighbours.
        #[test]
        fn card_bound_upper_bounds_tanimoto(a in arb_fp(), b in arb_fp()) {
            let bound = card_bound(a.magnitude(), b.magnitude());
            prop_assert!(bound >= a.tanimoto(&b) - 1e-9);
        }
    }
}
