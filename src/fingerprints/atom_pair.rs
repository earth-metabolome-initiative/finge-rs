use alloc::{vec, vec::Vec};

use crate::{
    bit_fingerprint::BitFingerprint,
    fingerprint::Fingerprint,
    traits::{AtomPairGraph, MolecularAtom, MolecularBond},
};

const COUNT_BOUNDS: [usize; 4] = [1, 2, 4, 8];

/// Bit-based 2D AtomPair fingerprint with RDKit-default count simulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomPairFingerprint {
    min_distance: u8,
    max_distance: u8,
    fp_size: usize,
    count_simulation: bool,
}

impl AtomPairFingerprint {
    /// Creates a new AtomPair fingerprint with RDKit-default distances.
    #[inline]
    #[must_use]
    pub const fn new(fp_size: usize) -> Self {
        Self {
            min_distance: 1,
            max_distance: 30,
            fp_size,
            count_simulation: true,
        }
    }

    /// Sets the minimum topological distance between paired atoms.
    #[inline]
    #[must_use]
    pub const fn with_min_distance(mut self, min_distance: u8) -> Self {
        self.min_distance = min_distance;
        self
    }

    /// Sets the maximum topological distance between paired atoms.
    #[inline]
    #[must_use]
    pub const fn with_max_distance(mut self, max_distance: u8) -> Self {
        self.max_distance = max_distance;
        self
    }

    /// Toggles RDKit-style count simulation.
    #[inline]
    #[must_use]
    pub const fn with_count_simulation(mut self, count_simulation: bool) -> Self {
        self.count_simulation = count_simulation;
        self
    }

    /// Returns the minimum topological distance.
    #[inline]
    #[must_use]
    pub const fn min_distance(self) -> u8 {
        self.min_distance
    }

    /// Returns the maximum topological distance.
    #[inline]
    #[must_use]
    pub const fn max_distance(self) -> u8 {
        self.max_distance
    }

    /// Returns the folded bit-vector length.
    #[inline]
    #[must_use]
    pub const fn fp_size(self) -> usize {
        self.fp_size
    }

    /// Returns whether count simulation is enabled.
    #[inline]
    #[must_use]
    pub const fn count_simulation(self) -> bool {
        self.count_simulation
    }
}

impl Default for AtomPairFingerprint {
    #[inline]
    fn default() -> Self {
        Self::new(2048)
    }
}

impl<G> Fingerprint<G> for AtomPairFingerprint
where
    G: AtomPairGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
{
    type Output = BitFingerprint;

    fn compute(&self, graph: &G) -> Self::Output {
        let mut fingerprint = BitFingerprint::zeros(self.fp_size);
        if self.fp_size == 0 || graph.atom_count() < 2 {
            return fingerprint;
        }

        let (adjacency, atom_codes) = graph.atom_pair_adjacency_and_codes();

        if self.count_simulation {
            let effective_size = self.fp_size / COUNT_BOUNDS.len();
            if effective_size == 0 {
                return fingerprint;
            }

            let mut counts = vec![0_u8; effective_size];
            if effective_size.is_power_of_two() {
                let mask = effective_size - 1;
                visit_atom_pairs(
                    &adjacency,
                    self.min_distance,
                    self.max_distance,
                    |left, right, distance| {
                        let bit_id =
                            hashed_atom_pair_bit(atom_codes[left], atom_codes[right], distance)
                                & mask;
                        counts[bit_id] = counts[bit_id].saturating_add(1).min(8);
                    },
                );
            } else {
                visit_atom_pairs(
                    &adjacency,
                    self.min_distance,
                    self.max_distance,
                    |left, right, distance| {
                        let bit_id =
                            hashed_atom_pair_bit(atom_codes[left], atom_codes[right], distance)
                                % effective_size;
                        counts[bit_id] = counts[bit_id].saturating_add(1).min(8);
                    },
                );
            }

            for (bit_id, count) in counts.into_iter().enumerate() {
                for (offset, bound) in COUNT_BOUNDS.into_iter().enumerate() {
                    if usize::from(count) >= bound {
                        fingerprint.set(bit_id * COUNT_BOUNDS.len() + offset);
                    }
                }
            }
        } else {
            if self.fp_size.is_power_of_two() {
                let mask = self.fp_size - 1;
                visit_atom_pairs(
                    &adjacency,
                    self.min_distance,
                    self.max_distance,
                    |left, right, distance| {
                        let bit_id =
                            hashed_atom_pair_bit(atom_codes[left], atom_codes[right], distance)
                                & mask;
                        fingerprint.set(bit_id);
                    },
                );
            } else {
                visit_atom_pairs(
                    &adjacency,
                    self.min_distance,
                    self.max_distance,
                    |left, right, distance| {
                        let bit_id =
                            hashed_atom_pair_bit(atom_codes[left], atom_codes[right], distance)
                                % self.fp_size;
                        fingerprint.set(bit_id);
                    },
                );
            }
        }

        fingerprint
    }
}

fn visit_atom_pairs<F>(adjacency: &[Vec<usize>], min_distance: u8, max_distance: u8, mut visit: F)
where
    F: FnMut(usize, usize, u8),
{
    let atom_count = adjacency.len();
    let mut distances = vec![u8::MAX; atom_count];
    let mut queue = Vec::with_capacity(atom_count);

    for start in 0..atom_count {
        distances.fill(u8::MAX);
        queue.clear();
        distances[start] = 0;
        queue.push(start);

        let mut head = 0;
        while head < queue.len() {
            let current = queue[head];
            head += 1;
            let next_distance = distances[current] + 1;

            if next_distance > max_distance {
                continue;
            }

            for &neighbor in &adjacency[current] {
                if distances[neighbor] != u8::MAX {
                    continue;
                }
                distances[neighbor] = next_distance;
                queue.push(neighbor);
            }
        }

        for (end, &distance) in distances.iter().enumerate().skip(start + 1) {
            if distance == u8::MAX {
                continue;
            }
            if distance < min_distance || distance > max_distance {
                continue;
            }
            visit(start, end, distance);
        }
    }
}

#[inline]
fn hashed_atom_pair_bit(left_code: u32, right_code: u32, distance: u8) -> usize {
    let (first, second) = if left_code <= right_code {
        (left_code, right_code)
    } else {
        (right_code, left_code)
    };

    let mut seed = 0_u32;
    hash_combine(&mut seed, first);
    hash_combine(&mut seed, u32::from(distance));
    hash_combine(&mut seed, second);
    seed as usize
}

#[inline]
fn hash_combine(seed: &mut u32, value: u32) {
    *seed ^= value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

#[cfg(test)]
mod tests {
    use alloc::{vec, vec::Vec};

    use smiles_parser::smiles::Smiles;

    use super::AtomPairFingerprint;
    use crate::{
        Fingerprint, smiles_support_impl::SmilesRdkitScratch,
        test_fixtures::rdkit_atom_pair_fixture,
    };

    fn observed_active_bits(smiles: &str, fingerprint: AtomPairFingerprint) -> Vec<usize> {
        let smiles: Smiles = smiles.parse().expect("fixture SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);
        fingerprint.compute(&graph).active_bits().collect()
    }

    #[test]
    fn rdkit_default_atom_pair_bit_fixtures_match() {
        for (smiles, expected_bits) in [
            ("CCC", vec![1148, 1404, 1405]),
            ("CCCC", vec![880, 1144, 1145, 1336, 1404, 1405]),
        ] {
            let observed = observed_active_bits(smiles, AtomPairFingerprint::default());
            assert_eq!(observed, expected_bits, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_atom_pair_bit_fixtures_without_count_simulation_match() {
        for (smiles, expected_bits) in [
            ("CCC", vec![1311, 1375]),
            ("CCCC", vec![1310, 1358, 1375, 1756]),
        ] {
            let observed = observed_active_bits(
                smiles,
                AtomPairFingerprint::default().with_count_simulation(false),
            );
            assert_eq!(observed, expected_bits, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_atom_pair_matrix_matches_reference_corpus() {
        let fixture = rdkit_atom_pair_fixture();
        assert_eq!(fixture.molecules.len(), 1_024);
        assert_eq!(fixture.cases.len(), 7);
        assert_eq!(
            fixture.source.dataset,
            "scikit-fingerprints HIV test corpus"
        );
        assert_eq!(
            fixture.source.selection,
            "1024 parseable SMILES fixture in repo order"
        );
        assert_eq!(fixture.source.generator, "RDKit AtomPairGenerator");
        assert!(!fixture.source.include_chirality);
        assert!(fixture.source.count_simulation);
        assert_eq!(fixture.source.min_distance, 1);
        assert_eq!(fixture.source.max_distance, 30);
        assert!(fixture.source.use_2d);
        assert_eq!(
            fixture
                .cases
                .iter()
                .map(|case| case.fp_size)
                .collect::<Vec<_>>(),
            [64, 128, 256, 512, 1024, 2048, 4096]
        );

        for case in &fixture.cases {
            for (smiles, expected_bits) in fixture.molecules.iter().zip(&case.active_bits) {
                let observed_bits =
                    observed_active_bits(smiles, AtomPairFingerprint::new(case.fp_size));
                assert_eq!(
                    observed_bits, *expected_bits,
                    "failed for fp_size={}, smiles={smiles}",
                    case.fp_size,
                );
            }
        }
    }
}
