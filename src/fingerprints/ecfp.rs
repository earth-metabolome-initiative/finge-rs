use alloc::{vec, vec::Vec};

use crate::{
    bit_fingerprint::BitFingerprint,
    fingerprint::Fingerprint,
    traits::{EcfpGraph, MolecularAtom, MolecularBond},
};

/// Bit-based extended-connectivity fingerprint (ECFP).
///
/// This implements the default Morgan/ECFP-style connectivity invariants and
/// produces a folded bit vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EcfpFingerprint {
    radius: u8,
    fp_size: usize,
    use_bond_types: bool,
    include_ring_membership: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct NeighborInfo {
    other: usize,
    edge_idx: usize,
    bond_invariant: u32,
}

impl EcfpFingerprint {
    /// Creates a new bit fingerprint with the requested radius and size.
    #[inline]
    #[must_use]
    pub const fn new(radius: u8, fp_size: usize) -> Self {
        Self {
            radius,
            fp_size,
            use_bond_types: true,
            include_ring_membership: true,
        }
    }

    /// Toggles bond-order-sensitive invariants.
    #[inline]
    #[must_use]
    pub const fn with_use_bond_types(mut self, use_bond_types: bool) -> Self {
        self.use_bond_types = use_bond_types;
        self
    }

    /// Toggles ring-membership participation in the initial atom invariants.
    #[inline]
    #[must_use]
    pub const fn with_include_ring_membership(mut self, include_ring_membership: bool) -> Self {
        self.include_ring_membership = include_ring_membership;
        self
    }

    /// Returns the ECFP radius.
    #[inline]
    #[must_use]
    pub const fn radius(self) -> u8 {
        self.radius
    }

    /// Returns the folded bit-vector length.
    #[inline]
    #[must_use]
    pub const fn fp_size(self) -> usize {
        self.fp_size
    }

    fn emit_hashes<G, F>(&self, graph: &G, mut emit_hash: F)
    where
        G: EcfpGraph<NodeId = usize>,
        G::NodeSymbol: MolecularAtom,
        G::Bond: MolecularBond<NodeId = usize>,
        F: FnMut(u32),
    {
        let atom_count = graph.atom_count();
        if atom_count == 0 {
            return;
        }

        let mut current_invariants = Vec::with_capacity(atom_count);

        for atom_id in 0..atom_count {
            let invariant = graph.ecfp_atom_invariant(atom_id, self.include_ring_membership);
            current_invariants.push(invariant);
            emit_hash(invariant);
        }

        let adjacency = adjacency(graph, self.use_bond_types);
        let mut next_layer_invariants = vec![0_u32; atom_count];
        let mut seen_neighborhoods = Vec::<Vec<usize>>::new();
        let mut current_atom_neighborhoods = vec![Vec::<usize>::new(); atom_count];
        let mut next_atom_neighborhoods = vec![Vec::<usize>::new(); atom_count];
        let mut dead_atoms = vec![false; atom_count];
        let mut ranked_atoms = Vec::<(usize, u32)>::with_capacity(atom_count);

        for layer in 0..self.radius {
            ranked_atoms.clear();

            for atom_id in 0..atom_count {
                if dead_atoms[atom_id] {
                    continue;
                }

                let neighbors = &adjacency[atom_id];
                if neighbors.is_empty() {
                    dead_atoms[atom_id] = true;
                    continue;
                }

                let mut neighborhood = core::mem::take(&mut next_atom_neighborhoods[atom_id]);
                neighborhood.clear();
                neighborhood.extend_from_slice(&current_atom_neighborhoods[atom_id]);
                let mut neighbor_pairs = Vec::with_capacity(neighbors.len());

                for neighbor in neighbors {
                    extend_atom_environment(
                        neighbor,
                        &current_atom_neighborhoods,
                        &mut neighborhood,
                        &current_invariants,
                        &mut neighbor_pairs,
                    );
                }

                neighbor_pairs.sort_unstable();
                let invariant = layered_invariant(
                    u32::from(layer),
                    current_invariants[atom_id],
                    &neighbor_pairs,
                );

                next_layer_invariants[atom_id] = invariant;
                next_atom_neighborhoods[atom_id] = neighborhood;
                ranked_atoms.push((atom_id, invariant));
            }

            ranked_atoms.sort_unstable_by(
                |(left_atom_id, left_invariant), (right_atom_id, right_invariant)| {
                    next_atom_neighborhoods[*left_atom_id]
                        .cmp(&next_atom_neighborhoods[*right_atom_id])
                        .then_with(|| left_invariant.cmp(right_invariant))
                        .then_with(|| left_atom_id.cmp(right_atom_id))
                },
            );
            for (atom_id, invariant) in &ranked_atoms {
                if mark_neighborhood_seen(
                    &mut seen_neighborhoods,
                    &next_atom_neighborhoods[*atom_id],
                ) {
                    emit_hash(*invariant);
                } else {
                    dead_atoms[*atom_id] = true;
                }
            }

            core::mem::swap(&mut current_invariants, &mut next_layer_invariants);
            next_layer_invariants.fill(0);
            core::mem::swap(
                &mut current_atom_neighborhoods,
                &mut next_atom_neighborhoods,
            );
        }
    }
}

impl Default for EcfpFingerprint {
    #[inline]
    fn default() -> Self {
        Self::new(2, 2048)
    }
}

impl<G> Fingerprint<G> for EcfpFingerprint
where
    G: EcfpGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
{
    type Output = BitFingerprint;

    fn compute(&self, graph: &G) -> Self::Output {
        let mut fingerprint = BitFingerprint::zeros(self.fp_size);
        if self.fp_size == 0 {
            return fingerprint;
        }

        if self.fp_size.is_power_of_two() {
            let mask = self.fp_size - 1;
            self.emit_hashes(graph, |hash| fingerprint.set(hash as usize & mask));
        } else {
            let fp_size = self.fp_size;
            self.emit_hashes(graph, |hash| fingerprint.set(hash as usize % fp_size));
        }

        fingerprint
    }
}

#[inline]
fn extend_atom_environment(
    neighbor: &NeighborInfo,
    atom_neighborhoods: &[Vec<usize>],
    neighborhood: &mut Vec<usize>,
    current_invariants: &[u32],
    neighbor_pairs: &mut Vec<(u32, u32)>,
) {
    insert_sorted_unique(neighborhood, neighbor.edge_idx);
    union_sorted_in_place(neighborhood, &atom_neighborhoods[neighbor.other]);
    neighbor_pairs.push((neighbor.bond_invariant, current_invariants[neighbor.other]));
}

#[inline]
fn layered_invariant(layer: u32, center_invariant: u32, neighbor_pairs: &[(u32, u32)]) -> u32 {
    let mut seed = layer;
    hash_combine(&mut seed, center_invariant);

    for &(bond_invariant, atom_invariant) in neighbor_pairs {
        hash_combine(&mut seed, pair_hash(bond_invariant, atom_invariant));
    }

    seed
}

#[inline]
fn pair_hash(left: u32, right: u32) -> u32 {
    let mut seed = 0_u32;
    hash_combine(&mut seed, left);
    hash_combine(&mut seed, right);
    seed
}

#[inline]
fn hash_combine(seed: &mut u32, value: u32) {
    *seed ^= value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

fn adjacency<G>(graph: &G, use_bond_types: bool) -> Vec<Vec<NeighborInfo>>
where
    G: EcfpGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
{
    let atom_count = graph.atom_count();
    let mut adjacency = vec![Vec::new(); atom_count];
    let mut edge_idx = 0;

    for node_id in 0..graph.atom_count() {
        for bond in graph.bonds(node_id) {
            let Some(other) = other_node(&bond, node_id) else {
                continue;
            };
            if other > node_id {
                let bond_invariant = graph.ecfp_bond_invariant(&bond, use_bond_types);
                adjacency[node_id].push(NeighborInfo {
                    other,
                    edge_idx,
                    bond_invariant,
                });
                adjacency[other].push(NeighborInfo {
                    other: node_id,
                    edge_idx,
                    bond_invariant,
                });
                edge_idx += 1;
            }
        }
    }

    adjacency
}

#[inline]
fn insert_sorted_unique(values: &mut Vec<usize>, value: usize) {
    match values.binary_search(&value) {
        Ok(_) => {}
        Err(index) => values.insert(index, value),
    }
}

#[inline]
fn union_sorted_in_place(target: &mut Vec<usize>, source: &[usize]) {
    if source.is_empty() {
        return;
    }

    let mut merged = Vec::with_capacity(target.len() + source.len());

    let mut left = 0;
    let mut right = 0;

    while left < target.len() && right < source.len() {
        match target[left].cmp(&source[right]) {
            core::cmp::Ordering::Less => {
                merged.push(target[left]);
                left += 1;
            }
            core::cmp::Ordering::Greater => {
                merged.push(source[right]);
                right += 1;
            }
            core::cmp::Ordering::Equal => {
                merged.push(target[left]);
                left += 1;
                right += 1;
            }
        }
    }

    merged.extend_from_slice(&target[left..]);
    merged.extend_from_slice(&source[right..]);
    *target = merged;
}

#[inline]
fn mark_neighborhood_seen(
    seen_neighborhoods: &mut Vec<Vec<usize>>,
    neighborhood: &[usize],
) -> bool {
    match seen_neighborhoods.binary_search_by(|candidate| candidate.as_slice().cmp(neighborhood)) {
        Ok(_) => false,
        Err(index) => {
            seen_neighborhoods.insert(index, neighborhood.to_vec());
            true
        }
    }
}

#[inline]
fn other_node<B>(bond: &B, atom_id: usize) -> Option<usize>
where
    B: MolecularBond<NodeId = usize>,
{
    if bond.source() == atom_id {
        Some(bond.target())
    } else if bond.target() == atom_id {
        Some(bond.source())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use alloc::{vec, vec::Vec};

    use smiles_parser::smiles::Smiles;

    use crate::smiles_support_impl::SmilesRdkitScratch;
    use crate::{Fingerprint, fingerprints::EcfpFingerprint, test_fixtures::rdkit_ecfp_fixture};

    fn observed_active_bits(smiles: &str, fingerprint: EcfpFingerprint) -> Vec<usize> {
        let smiles: Smiles = smiles.parse().expect("fixture SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);

        fingerprint.compute(&graph).active_bits().collect()
    }

    #[test]
    fn rdkit_ecfp4_bit_fixtures_match() {
        for (smiles, expected_bits) in [
            ("CCO", vec![80, 222, 294, 807, 1057, 1410]),
            ("CC=O", vec![308, 650, 694, 844, 1004, 1057]),
            ("C#N", vec![489, 915, 1384]),
            ("c1ccccc1", vec![389, 1088, 1873]),
            ("C1CCCCC1", vec![2, 926, 1028]),
            ("CC(C)O", vec![1, 227, 283, 709, 807, 1057]),
            ("[13CH3][NH3+]", vec![397, 1057, 1082]),
        ] {
            let observed_bits = observed_active_bits(smiles, EcfpFingerprint::default());
            assert_eq!(observed_bits, expected_bits, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_ecfp4_bit_fixtures_without_bond_types_match() {
        for (smiles, expected_bits) in [
            ("CC=O", vec![650, 694, 844, 1057, 1075, 1655]),
            ("C#N", vec![96, 915, 1384]),
        ] {
            let observed_bits = observed_active_bits(
                smiles,
                EcfpFingerprint::default().with_use_bond_types(false),
            );
            assert_eq!(observed_bits, expected_bits, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_ecfp_matrix_matches_reference_corpus() {
        let fixture = rdkit_ecfp_fixture();
        assert_eq!(fixture.molecules.len(), 100_000);
        assert_eq!(fixture.cases.len(), 42);
        assert_eq!(fixture.source.dataset, "PubChem CID-SMILES");
        assert_eq!(
            fixture.source.selection,
            "first 100000 unique SMILES parseable by smiles-parser and RDKit, in dataset order"
        );
        assert_eq!(fixture.source.generator, "RDKit MorganGenerator");
        assert!(!fixture.source.include_chirality);
        assert!(fixture.source.use_bond_types);
        assert!(fixture.source.include_ring_membership);
        assert_eq!(
            fixture
                .cases
                .iter()
                .map(|case| case.radius)
                .collect::<Vec<_>>(),
            [
                0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
                4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
            ]
        );
        assert_eq!(
            fixture
                .cases
                .iter()
                .map(|case| case.fp_size)
                .collect::<Vec<_>>(),
            [
                64, 128, 256, 512, 1024, 2048, 4096, 64, 128, 256, 512, 1024, 2048, 4096, 64, 128,
                256, 512, 1024, 2048, 4096, 64, 128, 256, 512, 1024, 2048, 4096, 64, 128, 256, 512,
                1024, 2048, 4096, 64, 128, 256, 512, 1024, 2048, 4096,
            ]
        );

        for case in &fixture.cases {
            for (smiles, expected_bits) in fixture.molecules.iter().zip(&case.active_bits) {
                let observed_bits =
                    observed_active_bits(smiles, EcfpFingerprint::new(case.radius, case.fp_size));

                assert_eq!(
                    observed_bits, *expected_bits,
                    "failed for radius={}, fp_size={}, smiles={smiles}",
                    case.radius, case.fp_size,
                );
            }
        }
    }

    #[test]
    fn preparing_the_same_smiles_twice_is_stable() {
        let smiles: Smiles = "c1ccncc1".parse().expect("fixture SMILES should parse");
        let fingerprint = EcfpFingerprint::default();
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);
        let reused = fingerprint.compute(&graph);

        let mut other_scratch = SmilesRdkitScratch::default();
        let other_graph = other_scratch.prepare(&smiles);
        let other = fingerprint.compute(&other_graph);

        assert_eq!(reused, other);
    }
}
