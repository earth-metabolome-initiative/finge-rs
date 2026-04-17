use alloc::{vec, vec::Vec};

/// Dense folded count fingerprint backed by `Vec<u32>`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CountFingerprint {
    counts: Vec<u32>,
}

impl CountFingerprint {
    /// Creates an all-zero count fingerprint with the requested length.
    #[inline]
    #[must_use]
    pub fn zeros(len: usize) -> Self {
        Self {
            counts: vec![0; len],
        }
    }

    /// Returns the raw count slice.
    #[inline]
    #[must_use]
    pub fn as_slice(&self) -> &[u32] {
        &self.counts
    }

    /// Increments one count when the index is in range.
    #[inline]
    pub fn increment(&mut self, index: usize) {
        if let Some(count) = self.counts.get_mut(index) {
            *count = count.saturating_add(1);
        }
    }

    /// Returns the count stored at `index`.
    #[inline]
    #[must_use]
    pub fn count(&self, index: usize) -> u32 {
        self.counts.get(index).copied().unwrap_or(0)
    }

    /// Returns the indices and counts for every nonzero entry.
    #[inline]
    pub fn active_counts(&self) -> impl Iterator<Item = (usize, u32)> + '_ {
        self.counts
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, count)| count != 0)
    }

    /// Returns the number of stored counts.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns whether the fingerprint stores no counts.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }
}

impl From<Vec<u32>> for CountFingerprint {
    #[inline]
    fn from(counts: Vec<u32>) -> Self {
        Self { counts }
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;

    use super::CountFingerprint;

    #[test]
    fn count_fingerprint_tracks_incremented_counts() {
        let mut fingerprint = CountFingerprint::zeros(8);
        fingerprint.increment(2);
        fingerprint.increment(2);
        fingerprint.increment(5);
        fingerprint.increment(99);

        assert_eq!(
            fingerprint.active_counts().collect::<alloc::vec::Vec<_>>(),
            vec![(2, 2), (5, 1)]
        );
        assert_eq!(fingerprint.count(2), 2);
        assert_eq!(fingerprint.count(3), 0);
        assert_eq!(fingerprint.len(), 8);
    }
}
