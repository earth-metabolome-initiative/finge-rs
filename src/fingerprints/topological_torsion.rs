use alloc::{collections::BTreeMap, vec, vec::Vec};

use crate::{
    bit_fingerprint::BitFingerprint,
    fingerprint::Fingerprint,
    traits::{MolecularAtom, MolecularBond, TopologicalTorsionGraph},
};

const COUNT_BOUNDS: [usize; 4] = [1, 2, 4, 8];
const MAX_TORSION_ATOM_COUNT: u8 = 7;

/// Bit-based RDKit-style Topological Torsion fingerprint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TopologicalTorsionFingerprint {
    torsion_atom_count: u8,
    fp_size: usize,
    count_simulation: bool,
    only_shortest_paths: bool,
}

impl TopologicalTorsionFingerprint {
    /// Creates a new Topological Torsion fingerprint with RDKit-default settings.
    #[inline]
    #[must_use]
    pub const fn new(fp_size: usize) -> Self {
        Self {
            torsion_atom_count: 4,
            fp_size,
            count_simulation: true,
            only_shortest_paths: false,
        }
    }

    /// Sets the number of atoms in each torsion path.
    #[inline]
    #[must_use]
    pub const fn with_torsion_atom_count(mut self, torsion_atom_count: u8) -> Self {
        self.torsion_atom_count = torsion_atom_count;
        self
    }

    /// Toggles RDKit-style count simulation.
    #[inline]
    #[must_use]
    pub const fn with_count_simulation(mut self, count_simulation: bool) -> Self {
        self.count_simulation = count_simulation;
        self
    }

    /// Toggles RDKit's shortest-path-only filtering.
    #[inline]
    #[must_use]
    pub const fn with_only_shortest_paths(mut self, only_shortest_paths: bool) -> Self {
        self.only_shortest_paths = only_shortest_paths;
        self
    }

    /// Returns the folded bit-vector length.
    #[inline]
    #[must_use]
    pub const fn fp_size(self) -> usize {
        self.fp_size
    }

    /// Returns the number of atoms in each torsion path.
    #[inline]
    #[must_use]
    pub const fn torsion_atom_count(self) -> u8 {
        self.torsion_atom_count
    }

    /// Returns whether count simulation is enabled.
    #[inline]
    #[must_use]
    pub const fn count_simulation(self) -> bool {
        self.count_simulation
    }

    /// Returns whether only shortest paths are included.
    #[inline]
    #[must_use]
    pub const fn only_shortest_paths(self) -> bool {
        self.only_shortest_paths
    }
}

impl Default for TopologicalTorsionFingerprint {
    #[inline]
    fn default() -> Self {
        Self::new(2048)
    }
}

impl<G> Fingerprint<G> for TopologicalTorsionFingerprint
where
    G: TopologicalTorsionGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
{
    type Output = BitFingerprint;

    fn compute(&self, graph: &G) -> Self::Output {
        let mut fingerprint = BitFingerprint::zeros(self.fp_size);
        if self.fp_size == 0 || self.torsion_atom_count < 2 || graph.atom_count() < 2 {
            return fingerprint;
        }
        assert!(
            self.torsion_atom_count <= MAX_TORSION_ATOM_COUNT,
            "TopologicalTorsionFingerprint supports torsion_atom_count <= 7 in non-chiral mode, matching RDKit's 64-bit sparse ids"
        );
        if graph.atom_count() < usize::from(self.torsion_atom_count) {
            return fingerprint;
        }

        let (adjacency, edge_count) = topological_torsion_edge_adjacency(graph);
        let endpoint_codes = (0..graph.atom_count())
            .map(|atom_id| graph.topological_torsion_atom_code(atom_id, 1))
            .collect::<Vec<_>>();
        let internal_codes = (0..graph.atom_count())
            .map(|atom_id| graph.topological_torsion_atom_code(atom_id, 2))
            .collect::<Vec<_>>();

        if self.count_simulation {
            let effective_size = self.fp_size / COUNT_BOUNDS.len();
            if effective_size == 0 {
                return fingerprint;
            }

            let mut counts = vec![0_u8; effective_size];
            if effective_size.is_power_of_two() {
                let mask = effective_size - 1;
                visit_topological_torsions(
                    &adjacency,
                    edge_count,
                    usize::from(self.torsion_atom_count),
                    self.only_shortest_paths,
                    |path| {
                        let bit_id =
                            hashed_topological_torsion_bit(path, &endpoint_codes, &internal_codes)
                                & mask;
                        counts[bit_id] = counts[bit_id].saturating_add(1).min(8);
                    },
                );
            } else {
                visit_topological_torsions(
                    &adjacency,
                    edge_count,
                    usize::from(self.torsion_atom_count),
                    self.only_shortest_paths,
                    |path| {
                        let bit_id =
                            hashed_topological_torsion_bit(path, &endpoint_codes, &internal_codes)
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
        } else if self.fp_size.is_power_of_two() {
            let mask = self.fp_size - 1;
            visit_topological_torsions(
                &adjacency,
                edge_count,
                usize::from(self.torsion_atom_count),
                self.only_shortest_paths,
                |path| {
                    let bit_id =
                        hashed_topological_torsion_bit(path, &endpoint_codes, &internal_codes)
                            & mask;
                    fingerprint.set(bit_id);
                },
            );
        } else {
            visit_topological_torsions(
                &adjacency,
                edge_count,
                usize::from(self.torsion_atom_count),
                self.only_shortest_paths,
                |path| {
                    let bit_id =
                        hashed_topological_torsion_bit(path, &endpoint_codes, &internal_codes)
                            % self.fp_size;
                    fingerprint.set(bit_id);
                },
            );
        }

        fingerprint
    }
}

fn topological_torsion_edge_adjacency<G>(graph: &G) -> (Vec<Vec<(usize, usize)>>, usize)
where
    G: TopologicalTorsionGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
{
    let atom_count = graph.atom_count();
    let mut adjacency = vec![Vec::new(); atom_count];
    let mut edge_ids = BTreeMap::new();
    let mut next_edge_id = 0_usize;

    for (atom_id, neighbors) in adjacency.iter_mut().enumerate() {
        for bond in graph.bonds(atom_id) {
            let (source, target) = (bond.source(), bond.target());
            let Some(neighbor) = (if source == atom_id {
                Some(target)
            } else if target == atom_id {
                Some(source)
            } else {
                None
            }) else {
                continue;
            };

            let key = if source <= target {
                (source, target)
            } else {
                (target, source)
            };
            let edge_id = *edge_ids.entry(key).or_insert_with(|| {
                let edge_id = next_edge_id;
                next_edge_id += 1;
                edge_id
            });
            neighbors.push((neighbor, edge_id));
        }
    }

    (adjacency, next_edge_id)
}

fn visit_topological_torsions<F>(
    adjacency: &[Vec<(usize, usize)>],
    edge_count: usize,
    torsion_atom_count: usize,
    only_shortest_paths: bool,
    mut visit: F,
) where
    F: FnMut(&[usize]),
{
    let atom_count = adjacency.len();
    let mut path = Vec::with_capacity(torsion_atom_count);
    let mut visited_edges = vec![false; edge_count];
    let mut shortest_distances = vec![usize::MAX; atom_count];

    for start in 0..atom_count {
        path.clear();
        visited_edges.fill(false);
        path.push(start);

        if only_shortest_paths {
            compute_shortest_distances(
                adjacency,
                start,
                torsion_atom_count - 1,
                &mut shortest_distances,
            );
        }

        extend_torsion_paths(
            adjacency,
            &mut path,
            &mut visited_edges,
            torsion_atom_count,
            if only_shortest_paths {
                Some(&shortest_distances)
            } else {
                None
            },
            &mut visit,
        );
    }
}

fn extend_torsion_paths<F>(
    adjacency: &[Vec<(usize, usize)>],
    path: &mut Vec<usize>,
    visited_edges: &mut [bool],
    torsion_atom_count: usize,
    shortest_distances: Option<&[usize]>,
    visit: &mut F,
) where
    F: FnMut(&[usize]),
{
    if path.len() == torsion_atom_count {
        if path.first() == path.last() {
            if !closed_cycle_path_is_canonical(path) {
                return;
            }
        } else if reverse_atom_path_is_preferred(path) {
            return;
        }
        if let Some(shortest_distances) = shortest_distances {
            let end = *path.last().expect("complete path should be non-empty");
            if shortest_distances[end] != torsion_atom_count - 1 {
                return;
            }
        }
        visit(path);
        return;
    }

    let current = *path.last().expect("path should be non-empty");
    for &(neighbor, edge_id) in &adjacency[current] {
        if visited_edges[edge_id] {
            continue;
        }
        visited_edges[edge_id] = true;
        path.push(neighbor);
        extend_torsion_paths(
            adjacency,
            path,
            visited_edges,
            torsion_atom_count,
            shortest_distances,
            visit,
        );
        path.pop();
        visited_edges[edge_id] = false;
    }
}

fn compute_shortest_distances(
    adjacency: &[Vec<(usize, usize)>],
    start: usize,
    max_distance: usize,
    distances: &mut [usize],
) {
    distances.fill(usize::MAX);
    distances[start] = 0;

    let mut queue = Vec::with_capacity(adjacency.len());
    queue.push(start);
    let mut head = 0;
    while head < queue.len() {
        let current = queue[head];
        head += 1;
        let next_distance = distances[current] + 1;
        if next_distance > max_distance {
            continue;
        }
        for &(neighbor, _) in &adjacency[current] {
            if distances[neighbor] != usize::MAX {
                continue;
            }
            distances[neighbor] = next_distance;
            queue.push(neighbor);
        }
    }
}

#[inline]
fn reverse_atom_path_is_preferred(path: &[usize]) -> bool {
    for index in 0..path.len() {
        let left = path[index];
        let right = path[path.len() - 1 - index];
        if left < right {
            return false;
        }
        if left > right {
            return true;
        }
    }
    false
}

fn closed_cycle_path_is_canonical(path: &[usize]) -> bool {
    debug_assert!(path.len() >= 2);
    debug_assert_eq!(path.first(), path.last());

    let cycle = &path[..path.len() - 1];
    let cycle_len = cycle.len();
    let mut best = Vec::with_capacity(path.len());
    best.extend_from_slice(path);

    for start in 0..cycle_len {
        let mut candidate = Vec::with_capacity(path.len());
        for offset in 0..cycle_len {
            candidate.push(cycle[(start + offset) % cycle_len]);
        }
        candidate.push(candidate[0]);
        if candidate < best {
            best = candidate;
        }
    }

    for start in 0..cycle_len {
        let mut candidate = Vec::with_capacity(path.len());
        for offset in 0..cycle_len {
            let index = (start + cycle_len - offset) % cycle_len;
            candidate.push(cycle[index]);
        }
        candidate.push(candidate[0]);
        if candidate < best {
            best = candidate;
        }
    }

    path == best.as_slice()
}

#[inline]
fn hashed_topological_torsion_bit(
    path: &[usize],
    endpoint_codes: &[u32],
    internal_codes: &[u32],
) -> usize {
    let forward_is_canonical =
        torsion_codes_are_lexicographically_leq(path, endpoint_codes, internal_codes);
    let mut seed = 0_u32;

    for out_index in 0..path.len() {
        let atom_index = if forward_is_canonical {
            path[out_index]
        } else {
            path[path.len() - 1 - out_index]
        };
        let atom_code = if out_index == 0 || out_index + 1 == path.len() {
            endpoint_codes[atom_index]
        } else {
            internal_codes[atom_index]
        };
        hash_combine(&mut seed, atom_code);
    }

    seed as usize
}

#[inline]
fn hash_combine(seed: &mut u32, value: u32) {
    *seed ^= value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

#[inline]
fn torsion_codes_are_lexicographically_leq(
    path: &[usize],
    endpoint_codes: &[u32],
    internal_codes: &[u32],
) -> bool {
    for index in 0..path.len() {
        let left = torsion_atom_code_for_position(path, index, endpoint_codes, internal_codes);
        let right = torsion_atom_code_for_position(
            path,
            path.len() - 1 - index,
            endpoint_codes,
            internal_codes,
        );
        if left < right {
            return true;
        }
        if left > right {
            return false;
        }
    }
    true
}

#[inline]
fn torsion_atom_code_for_position(
    path: &[usize],
    index: usize,
    endpoint_codes: &[u32],
    internal_codes: &[u32],
) -> u32 {
    let atom_id = path[index];
    if index == 0 || index + 1 == path.len() {
        endpoint_codes[atom_id]
    } else {
        internal_codes[atom_id]
    }
}

#[cfg(test)]
mod tests {
    use alloc::{vec, vec::Vec};

    use smiles_parser::smiles::Smiles;

    use super::TopologicalTorsionFingerprint;
    use crate::{
        Fingerprint, smiles_support_impl::SmilesRdkitScratch,
        test_fixtures::rdkit_topological_torsion_fixture,
    };

    fn observed_active_bits(
        smiles: &str,
        fingerprint: TopologicalTorsionFingerprint,
    ) -> Vec<usize> {
        let smiles: Smiles = smiles.parse().expect("fixture SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);
        fingerprint.compute(&graph).active_bits().collect()
    }

    #[test]
    fn rdkit_default_topological_torsion_bit_fixtures_match() {
        for (smiles, expected_bits) in [("CCCCC", vec![0, 1]), ("C1CCC1", vec![284, 285, 286])] {
            let observed = observed_active_bits(smiles, TopologicalTorsionFingerprint::default());
            assert_eq!(observed, expected_bits, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_topological_torsion_bit_fixtures_without_count_simulation_match() {
        for (smiles, expected_bits) in [("CCCCC", vec![0]), ("C1CCC1", vec![71])] {
            let observed = observed_active_bits(
                smiles,
                TopologicalTorsionFingerprint::default().with_count_simulation(false),
            );
            assert_eq!(observed, expected_bits, "failed for {smiles}");
        }
    }

    #[test]
    fn shortest_path_filter_matches_rdkit_ring_behavior() {
        let observed = observed_active_bits(
            "C1CCC1",
            TopologicalTorsionFingerprint::default()
                .with_count_simulation(false)
                .with_only_shortest_paths(true),
        );
        assert!(observed.is_empty());
    }

    #[test]
    fn rdkit_topological_torsion_matrix_matches_reference_corpus() {
        let fixture = rdkit_topological_torsion_fixture();
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
        assert_eq!(
            fixture.source.generator,
            "RDKit TopologicalTorsionGenerator"
        );
        assert!(!fixture.source.include_chirality);
        assert!(fixture.source.count_simulation);
        assert!(!fixture.source.only_shortest_paths);
        assert_eq!(fixture.source.torsion_atom_count, 4);
        assert_eq!(
            fixture
                .cases
                .iter()
                .map(|case| case.fp_size)
                .collect::<Vec<_>>(),
            [64, 128, 256, 512, 1024, 2048, 4096]
        );

        for case in &fixture.cases {
            let fingerprint = TopologicalTorsionFingerprint::new(case.fp_size);
            for (index, (smiles, expected_bits)) in
                fixture.molecules.iter().zip(&case.active_bits).enumerate()
            {
                let observed = observed_active_bits(smiles, fingerprint);
                assert_eq!(
                    observed, *expected_bits,
                    "failed for {} at molecule index {}",
                    smiles, index
                );
            }
        }
    }
}
