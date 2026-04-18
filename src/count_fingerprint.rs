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

/// Exact-radius layered count fingerprint backed by one [`CountFingerprint`]
/// per Morgan layer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LayeredCountFingerprint {
    layers: Vec<CountFingerprint>,
}

impl LayeredCountFingerprint {
    /// Creates an all-zero layered count fingerprint.
    #[inline]
    #[must_use]
    pub fn zeros(layer_count: usize, len: usize) -> Self {
        Self {
            layers: (0..layer_count)
                .map(|_| CountFingerprint::zeros(len))
                .collect(),
        }
    }

    /// Returns the per-layer fingerprints.
    #[inline]
    #[must_use]
    pub fn layers(&self) -> &[CountFingerprint] {
        &self.layers
    }

    /// Returns one layer by exact radius.
    #[inline]
    #[must_use]
    pub fn layer(&self, radius: usize) -> Option<&CountFingerprint> {
        self.layers.get(radius)
    }

    /// Returns the radius-0 layer.
    #[inline]
    #[must_use]
    pub fn formula(&self) -> &CountFingerprint {
        &self.layers[0]
    }

    /// Returns the number of exact-radius layers.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.layers.len()
    }

    /// Returns whether no layers are stored.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.layers.is_empty()
    }

    /// Increments one bit count in one exact-radius layer.
    #[inline]
    pub fn increment(&mut self, radius: usize, index: usize) {
        if let Some(layer) = self.layers.get_mut(radius) {
            layer.increment(index);
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;

    use super::{CountFingerprint, LayeredCountFingerprint};

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

    #[test]
    fn layered_count_fingerprint_tracks_per_radius_counts() {
        let mut fingerprint = LayeredCountFingerprint::zeros(3, 8);
        fingerprint.increment(0, 1);
        fingerprint.increment(2, 5);
        fingerprint.increment(2, 5);

        assert_eq!(fingerprint.len(), 3);
        assert_eq!(
            fingerprint
                .formula()
                .active_counts()
                .collect::<alloc::vec::Vec<_>>(),
            vec![(1, 1)]
        );
        assert_eq!(
            fingerprint
                .layer(2)
                .expect("layer 2 should exist")
                .active_counts()
                .collect::<alloc::vec::Vec<_>>(),
            vec![(5, 2)]
        );
    }
}
