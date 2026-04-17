use bitvec::prelude::{BitSlice, BitVec, Lsb0};

/// Dense bit fingerprint backed by `bitvec`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BitFingerprint {
    bits: BitVec<usize, Lsb0>,
}

impl BitFingerprint {
    /// Creates an all-zero fingerprint with the requested bit length.
    #[inline]
    #[must_use]
    pub fn zeros(len: usize) -> Self {
        Self {
            bits: BitVec::repeat(false, len),
        }
    }

    /// Returns the underlying bit slice.
    #[inline]
    #[must_use]
    pub fn as_bitslice(&self) -> &BitSlice<usize, Lsb0> {
        self.bits.as_bitslice()
    }

    /// Sets a bit to `true` when the index is in range.
    #[inline]
    pub fn set(&mut self, index: usize) {
        if index < self.bits.len() {
            self.bits.set(index, true);
        }
    }

    /// Returns whether the bit at `index` is set.
    #[inline]
    #[must_use]
    pub fn contains(&self, index: usize) -> bool {
        self.bits.get(index).is_some_and(|bit| *bit)
    }

    /// Returns the indices of every set bit.
    #[inline]
    pub fn active_bits(&self) -> impl Iterator<Item = usize> + '_ {
        self.bits.iter_ones()
    }

    /// Returns the number of bits in the fingerprint.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.bits.len()
    }

    /// Returns whether the fingerprint stores no bits.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bits.is_empty()
    }
}

impl From<BitVec<usize, Lsb0>> for BitFingerprint {
    #[inline]
    fn from(bits: BitVec<usize, Lsb0>) -> Self {
        Self { bits }
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;

    use super::BitFingerprint;

    #[test]
    fn bit_fingerprint_tracks_set_bits() {
        let mut fingerprint = BitFingerprint::zeros(16);
        fingerprint.set(1);
        fingerprint.set(7);
        fingerprint.set(99);

        assert_eq!(
            fingerprint.active_bits().collect::<alloc::vec::Vec<_>>(),
            vec![1, 7]
        );
        assert!(fingerprint.contains(1));
        assert!(!fingerprint.contains(2));
        assert_eq!(fingerprint.len(), 16);
    }
}
