use alloc::{vec, vec::Vec};

use crate::{
    bit_fingerprint::BitFingerprint,
    fingerprint::Fingerprint,
    traits::{MolecularAtom, MolecularBond, RdkFingerprintBond, RdkFingerprintGraph},
};

/// Bit-based RDKit-style topological path fingerprint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RdkFingerprint {
    min_path: u8,
    max_path: u8,
    fp_size: usize,
    n_bits_per_hash: u8,
    use_hs: bool,
    tgt_density: f64,
    min_size: usize,
    branched_paths: bool,
    use_bond_order: bool,
}

impl RdkFingerprint {
    /// Creates a new RDKFingerprint with RDKit-default settings.
    #[inline]
    #[must_use]
    pub const fn new(fp_size: usize) -> Self {
        Self {
            min_path: 1,
            max_path: 7,
            fp_size,
            n_bits_per_hash: 2,
            use_hs: true,
            tgt_density: 0.0,
            min_size: 128,
            branched_paths: true,
            use_bond_order: true,
        }
    }

    #[inline]
    #[must_use]
    pub const fn with_min_path(mut self, min_path: u8) -> Self {
        self.min_path = min_path;
        self
    }

    #[inline]
    #[must_use]
    pub const fn with_max_path(mut self, max_path: u8) -> Self {
        self.max_path = max_path;
        self
    }

    #[inline]
    #[must_use]
    pub const fn with_n_bits_per_hash(mut self, n_bits_per_hash: u8) -> Self {
        self.n_bits_per_hash = n_bits_per_hash;
        self
    }

    #[inline]
    #[must_use]
    pub const fn with_use_hs(mut self, use_hs: bool) -> Self {
        self.use_hs = use_hs;
        self
    }

    #[inline]
    #[must_use]
    pub const fn with_min_size(mut self, min_size: usize) -> Self {
        self.min_size = min_size;
        self
    }

    #[inline]
    #[must_use]
    pub const fn with_branched_paths(mut self, branched_paths: bool) -> Self {
        self.branched_paths = branched_paths;
        self
    }

    #[inline]
    #[must_use]
    pub const fn with_use_bond_order(mut self, use_bond_order: bool) -> Self {
        self.use_bond_order = use_bond_order;
        self
    }

    #[inline]
    #[must_use]
    pub const fn with_tgt_density(mut self, tgt_density: f64) -> Self {
        self.tgt_density = tgt_density;
        self
    }

    #[inline]
    #[must_use]
    pub const fn min_path(self) -> u8 {
        self.min_path
    }

    #[inline]
    #[must_use]
    pub const fn max_path(self) -> u8 {
        self.max_path
    }

    #[inline]
    #[must_use]
    pub const fn fp_size(self) -> usize {
        self.fp_size
    }

    #[inline]
    #[must_use]
    pub const fn n_bits_per_hash(self) -> u8 {
        self.n_bits_per_hash
    }

    #[inline]
    #[must_use]
    pub const fn use_hs(self) -> bool {
        self.use_hs
    }

    #[inline]
    #[must_use]
    pub const fn tgt_density(self) -> f64 {
        self.tgt_density
    }

    #[inline]
    #[must_use]
    pub const fn min_size(self) -> usize {
        self.min_size
    }

    #[inline]
    #[must_use]
    pub const fn branched_paths(self) -> bool {
        self.branched_paths
    }

    #[inline]
    #[must_use]
    pub const fn use_bond_order(self) -> bool {
        self.use_bond_order
    }
}

impl Default for RdkFingerprint {
    #[inline]
    fn default() -> Self {
        Self::new(2048)
    }
}

impl<G> Fingerprint<G> for RdkFingerprint
where
    G: RdkFingerprintGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
{
    type Output = BitFingerprint;

    fn compute(&self, graph: &G) -> Self::Output {
        if self.fp_size == 0
            || self.min_path == 0
            || self.max_path < self.min_path
            || self.n_bits_per_hash == 0
        {
            return BitFingerprint::zeros(self.fp_size);
        }

        let (atom_invariants, bonds, atom_bonds) =
            graph.rdk_fingerprint_graph(self.use_hs, self.use_bond_order);
        if bonds.is_empty() {
            return BitFingerprint::zeros(self.fp_size);
        }

        let mut fingerprint = BitFingerprint::zeros(self.fp_size);
        if self.branched_paths {
            visit_branched_subgraphs(
                &bonds,
                &atom_bonds,
                usize::from(self.min_path),
                usize::from(self.max_path),
                |bond_path| {
                    let seed = hashed_rdk_feature(bond_path, &bonds, &atom_invariants);
                    set_rdk_feature_bits(
                        &mut fingerprint,
                        seed,
                        self.fp_size,
                        self.n_bits_per_hash,
                    );
                },
            );
        } else {
            visit_linear_paths(
                &bonds,
                &atom_bonds,
                usize::from(self.min_path),
                usize::from(self.max_path),
                |bond_path| {
                    let seed = hashed_rdk_feature(bond_path, &bonds, &atom_invariants);
                    set_rdk_feature_bits(
                        &mut fingerprint,
                        seed,
                        self.fp_size,
                        self.n_bits_per_hash,
                    );
                },
            );
        }

        fold_to_density(fingerprint, self.tgt_density, self.min_size)
    }
}

fn visit_branched_subgraphs<F>(
    bonds: &[RdkFingerprintBond],
    atom_bonds: &[Vec<usize>],
    min_path: usize,
    max_path: usize,
    visit: F,
) where
    F: FnMut(&[usize]),
{
    BranchedWalker {
        bonds,
        atom_bonds,
        min_path,
        max_path,
        visit,
    }
    .run();
}

fn visit_linear_paths<F>(
    bonds: &[RdkFingerprintBond],
    atom_bonds: &[Vec<usize>],
    min_path: usize,
    max_path: usize,
    visit: F,
) where
    F: FnMut(&[usize]),
{
    LinearWalker {
        bonds,
        atom_bonds,
        min_path,
        max_path,
        visit,
    }
    .run();
}

struct BranchedWalker<'a, F> {
    bonds: &'a [RdkFingerprintBond],
    atom_bonds: &'a [Vec<usize>],
    min_path: usize,
    max_path: usize,
    visit: F,
}

impl<F> BranchedWalker<'_, F>
where
    F: FnMut(&[usize]),
{
    fn run(&mut self) {
        let mut global_forbidden = vec![false; self.bonds.len()];

        for start_bond in 0..self.bonds.len() {
            if global_forbidden[start_bond] {
                continue;
            }
            global_forbidden[start_bond] = true;

            let mut path = vec![start_bond];
            let mut candidates = bond_neighbors(start_bond, self.bonds, self.atom_bonds);
            candidates.sort_unstable();
            candidates.dedup();
            self.recurse(&mut path, &mut candidates, global_forbidden.clone());
        }
    }

    fn recurse(
        &mut self,
        path: &mut Vec<usize>,
        candidates: &mut Vec<usize>,
        forbidden: Vec<bool>,
    ) {
        let path_len = path.len();
        if path_len >= self.min_path && path_len <= self.max_path {
            (self.visit)(path);
        }
        if path_len >= self.max_path {
            return;
        }

        while let Some(next_bond) = candidates.pop() {
            if forbidden[next_bond] {
                continue;
            }

            let mut next_forbidden = forbidden.clone();
            next_forbidden[next_bond] = true;

            let mut next_candidates = candidates.clone();
            next_candidates.extend(
                bond_neighbors(next_bond, self.bonds, self.atom_bonds)
                    .into_iter()
                    .filter(|&bond_id| !next_forbidden[bond_id]),
            );
            next_candidates.sort_unstable();
            next_candidates.dedup();

            path.push(next_bond);
            self.recurse(path, &mut next_candidates, next_forbidden);
            path.pop();
        }
    }
}

struct LinearWalker<'a, F> {
    bonds: &'a [RdkFingerprintBond],
    atom_bonds: &'a [Vec<usize>],
    min_path: usize,
    max_path: usize,
    visit: F,
}

impl<F> LinearWalker<'_, F>
where
    F: FnMut(&[usize]),
{
    fn run(&mut self) {
        let atom_count = self.atom_bonds.len();
        let mut path_atoms = Vec::with_capacity(self.max_path + 1);
        let mut path_bonds = Vec::with_capacity(self.max_path);
        let mut visited_atoms = vec![false; atom_count];

        for start_atom in 0..atom_count {
            visited_atoms[start_atom] = true;
            path_atoms.push(start_atom);
            self.recurse(
                start_atom,
                &mut path_atoms,
                &mut path_bonds,
                &mut visited_atoms,
            );
            path_atoms.pop();
            visited_atoms[start_atom] = false;
        }
    }

    fn recurse(
        &mut self,
        current_atom: usize,
        path_atoms: &mut Vec<usize>,
        path_bonds: &mut Vec<usize>,
        visited_atoms: &mut [bool],
    ) {
        let path_len = path_bonds.len();
        if path_len >= self.min_path && path_len <= self.max_path {
            (self.visit)(path_bonds);
        }
        if path_len >= self.max_path {
            return;
        }

        for &bond_id in &self.atom_bonds[current_atom] {
            let bond = self.bonds[bond_id];
            let next_atom = if bond.source == current_atom {
                bond.target
            } else {
                bond.source
            };

            let allow_ring_closure = self.max_path > 2
                && path_len + 1 == self.max_path
                && path_atoms.len() > 2
                && next_atom == path_atoms[0]
                && path_atoms[path_atoms.len() - 2] != next_atom;

            if visited_atoms[next_atom] && !allow_ring_closure {
                continue;
            }

            path_bonds.push(bond_id);
            if allow_ring_closure {
                path_atoms.push(next_atom);
                (self.visit)(path_bonds);
                path_atoms.pop();
            } else {
                visited_atoms[next_atom] = true;
                path_atoms.push(next_atom);
                self.recurse(next_atom, path_atoms, path_bonds, visited_atoms);
                path_atoms.pop();
                visited_atoms[next_atom] = false;
            }
            path_bonds.pop();
        }
    }
}

fn bond_neighbors(
    bond_id: usize,
    bonds: &[RdkFingerprintBond],
    atom_bonds: &[Vec<usize>],
) -> Vec<usize> {
    let bond = bonds[bond_id];
    let mut neighbors = atom_bonds[bond.source].clone();
    neighbors.extend_from_slice(&atom_bonds[bond.target]);
    neighbors.retain(|&other_bond| other_bond != bond_id);
    neighbors
}

fn hashed_rdk_feature(
    bond_path: &[usize],
    bonds: &[RdkFingerprintBond],
    atom_invariants: &[u32],
) -> u32 {
    let atom_count = atom_invariants.len();
    let mut atom_degrees = vec![0_u32; atom_count];
    let mut atoms_in_path_count = 0_u32;
    for &bond_id in bond_path {
        let bond = bonds[bond_id];
        if atom_degrees[bond.source] == 0 {
            atoms_in_path_count += 1;
        }
        if atom_degrees[bond.target] == 0 {
            atoms_in_path_count += 1;
        }
        atom_degrees[bond.source] += 1;
        atom_degrees[bond.target] += 1;
    }

    let mut bond_hashes = Vec::with_capacity(bond_path.len());
    for &bond_id in bond_path {
        let bond = bonds[bond_id];
        let bond_neighbors = bond_path
            .iter()
            .copied()
            .filter(|&other_id| other_id != bond_id && bonds_share_atom(bond, bonds[other_id]))
            .count() as u32;

        let mut a1_hash = atom_invariants[bond.source];
        let mut a2_hash = atom_invariants[bond.target];
        let mut deg1 = atom_degrees[bond.source];
        let mut deg2 = atom_degrees[bond.target];
        if a1_hash < a2_hash {
            core::mem::swap(&mut a1_hash, &mut a2_hash);
            core::mem::swap(&mut deg1, &mut deg2);
        } else if a1_hash == a2_hash && deg1 < deg2 {
            core::mem::swap(&mut deg1, &mut deg2);
        }

        let mut bond_hash = bond_neighbors;
        hash_combine(&mut bond_hash, bond.bond_code);
        hash_combine(&mut bond_hash, a1_hash);
        hash_combine(&mut bond_hash, deg1);
        hash_combine(&mut bond_hash, a2_hash);
        hash_combine(&mut bond_hash, deg2);
        bond_hashes.push(bond_hash);
    }

    if bond_hashes.len() == 1 {
        return bond_hashes[0];
    }

    bond_hashes.sort_unstable();
    bond_hashes.push(atoms_in_path_count);
    hash_range(&bond_hashes)
}

#[inline]
fn bonds_share_atom(left: RdkFingerprintBond, right: RdkFingerprintBond) -> bool {
    left.source == right.source
        || left.source == right.target
        || left.target == right.source
        || left.target == right.target
}

#[inline]
fn hash_range(values: &[u32]) -> u32 {
    let mut seed = 0_u32;
    for &value in values {
        hash_combine(&mut seed, value);
    }
    seed
}

#[inline]
fn hash_combine(seed: &mut u32, value: u32) {
    *seed ^= value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

fn set_rdk_feature_bits(
    fingerprint: &mut BitFingerprint,
    seed: u32,
    fp_size: usize,
    n_bits_per_hash: u8,
) {
    let mut bit_id = seed as usize;
    if fp_size != 0 {
        bit_id %= fp_size;
    }
    fingerprint.set(bit_id);

    if n_bits_per_hash <= 1 {
        return;
    }

    let mut generator = RdkitMersenneTwister::seeded(seed);
    for _ in 1..n_bits_per_hash {
        let bit_id = (generator.next_u32() >> 1) as usize % fp_size;
        fingerprint.set(bit_id);
    }
}

fn fold_to_density(
    mut fingerprint: BitFingerprint,
    tgt_density: f64,
    min_size: usize,
) -> BitFingerprint {
    if tgt_density <= 0.0 {
        return fingerprint;
    }

    while fingerprint.len() >= min_size.saturating_mul(2) {
        let density = fingerprint.active_bits().count() as f64 / fingerprint.len() as f64;
        if density >= tgt_density {
            break;
        }

        let next_len = fingerprint.len() / 2;
        let mut folded = BitFingerprint::zeros(next_len);
        for bit in fingerprint.active_bits() {
            folded.set(bit % next_len);
        }
        fingerprint = folded;
    }

    fingerprint
}

#[derive(Clone, Copy)]
struct RdkitMersenneTwister {
    state: [u32; 4],
    index: usize,
}

impl RdkitMersenneTwister {
    const WORD_SIZE: u32 = 32;
    const UPPER_MASK: u32 = 0x8000_0000;
    const LOWER_MASK: u32 = 0x7fff_ffff;
    const MATRIX_A: u32 = 0x9908_b0df;
    const INIT_MULTIPLIER: u32 = 1_812_433_253;

    #[inline]
    fn seeded(seed: u32) -> Self {
        let mut state = [0_u32; 4];
        state[0] = seed;
        for i in 1..4 {
            state[i] = Self::INIT_MULTIPLIER
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }

        let mut generator = Self { state, index: 4 };
        generator.normalize_state();
        generator
    }

    #[inline]
    fn next_u32(&mut self) -> u32 {
        if self.index == 4 {
            self.twist();
        }

        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    fn normalize_state(&mut self) {
        let mut y0 = self.state[1] ^ self.state[3];
        if (y0 & (1_u32 << (Self::WORD_SIZE - 1))) != 0 {
            y0 = ((y0 ^ Self::MATRIX_A) << 1) | 1;
        } else {
            y0 <<= 1;
        }
        self.state[0] = (self.state[0] & Self::UPPER_MASK) | (y0 & Self::LOWER_MASK);

        if self.state.iter().all(|&value| value == 0) {
            self.state[0] = 1_u32 << (Self::WORD_SIZE - 1);
        }
    }

    fn twist(&mut self) {
        for i in 0..4 {
            let y =
                (self.state[i] & Self::UPPER_MASK) | (self.state[(i + 1) % 4] & Self::LOWER_MASK);
            let mut next = self.state[(i + 2) % 4] ^ (y >> 1);
            if (y & 1) != 0 {
                next ^= Self::MATRIX_A;
            }
            self.state[i] = next;
        }
        self.index = 0;
    }
}

#[cfg(test)]
mod tests {
    use alloc::{vec, vec::Vec};

    use elements_rs::Element;
    use smiles_parser::{atom::Atom, smiles::Smiles};

    use super::{RdkFingerprint, hashed_rdk_feature, set_rdk_feature_bits};
    use crate::{
        BitFingerprint, Fingerprint,
        smiles_support_impl::SmilesRdkitScratch,
        test_fixtures::rdkit_rdk_fixture,
        traits::{MolecularGraph, RdkFingerprintBond, RdkFingerprintGraph},
    };

    fn observed_active_bits(smiles: &str, fingerprint: RdkFingerprint) -> Vec<usize> {
        let smiles: Smiles = smiles.parse().expect("fixture SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);
        fingerprint.compute(&graph).active_bits().collect()
    }

    #[test]
    fn rdkit_default_rdk_bit_fixtures_match() {
        for (smiles, expected_bits) in [
            ("CCC", vec![709, 1308, 1772, 1813]),
            ("CCCC", vec![709, 875, 1308, 1772, 1813, 1927]),
            ("CCO", vec![562, 1183, 1308, 1339, 1728, 1772]),
        ] {
            let observed = observed_active_bits(smiles, RdkFingerprint::default());
            assert_eq!(observed, expected_bits, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_rdk_bit_fixtures_with_one_bit_per_hash_match() {
        for (smiles, expected_bits) in [
            ("CCC", vec![1308, 1813]),
            ("CCCC", vec![875, 1308, 1813]),
            ("CCO", vec![1308, 1339, 1728]),
        ] {
            let observed =
                observed_active_bits(smiles, RdkFingerprint::default().with_n_bits_per_hash(1));
            assert_eq!(observed, expected_bits, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_rdk_collapses_plain_explicit_hydrogens_on_heavy_atoms() {
        for (smiles, expected_bits) in [
            ("[H]C", vec![]),
            ("[2H]C", vec![324, 956]),
            ("[H]C([3H])", vec![324, 956]),
            ("[H]O[H]", vec![]),
            ("O=[ClH]", vec![148, 480]),
            ("[O-]Cl=O.[Na+]", vec![148, 302, 480, 543]),
            (
                "O=I(O)(O)(O)(O)O",
                vec![
                    504, 565, 809, 819, 1122, 1166, 1409, 1639, 1659, 1685, 1722, 1820,
                ],
            ),
        ] {
            let observed = observed_active_bits(smiles, RdkFingerprint::default());
            assert_eq!(observed, expected_bits, "failed for {smiles}");
        }
    }

    const FIRST_100K_COUNTEREXAMPLE: &str = "CN(C)C=CC=CC=[N+](C)C.[O-]Cl(=O)(=O)=O";
    const SECOND_500K_COUNTEREXAMPLE: &str = "O=C1C2=CC=CC=C2I(O)(O1)=O";
    const THIRD_500K_COUNTEREXAMPLE: &str = "O=C1CN2C(C3=CC=CC=C3I2(O1)=O)=O";
    const FIRST_1M_COUNTEREXAMPLE: &str = "CCCC[C@@](C)(C(=O)O)N=P(=O)OC1=C(C=CC2=CC=CC=C21)OC[C@@H]3[C@H]([C@@]([C@@H](O3)N4C5=C(C(=NC(=N5)N)OC)N=C4I)(C)O)O";

    #[test]
    fn rdk_first_100k_counterexample_matches_rdkit_default() {
        let observed = observed_active_bits(FIRST_100K_COUNTEREXAMPLE, RdkFingerprint::default());
        let expected = vec![
            14, 22, 80, 107, 118, 122, 148, 163, 166, 173, 194, 238, 244, 283, 302, 308, 313, 315,
            322, 352, 412, 427, 429, 433, 452, 456, 480, 543, 550, 592, 609, 621, 681, 694, 700,
            725, 750, 884, 894, 904, 935, 955, 956, 999, 1004, 1057, 1060, 1074, 1079, 1105, 1110,
            1127, 1142, 1147, 1194, 1220, 1241, 1271, 1284, 1294, 1308, 1325, 1339, 1425, 1444,
            1465, 1531, 1553, 1574, 1602, 1606, 1661, 1740, 1772, 1780, 1835, 1851, 1890, 1930,
            1950, 2004, 2019, 2047,
        ];
        assert_eq!(observed, expected);
    }

    #[test]
    fn rdk_first_100k_counterexample_matches_rdkit_with_one_bit_per_hash() {
        let observed = observed_active_bits(
            FIRST_100K_COUNTEREXAMPLE,
            RdkFingerprint::default().with_n_bits_per_hash(1),
        );
        let expected = vec![
            14, 22, 173, 194, 308, 313, 315, 322, 352, 412, 429, 452, 480, 543, 609, 681, 694, 725,
            894, 904, 935, 956, 999, 1057, 1105, 1127, 1142, 1147, 1220, 1241, 1308, 1325, 1425,
            1444, 1661, 1780, 1890, 1930, 2004, 2019,
        ];
        assert_eq!(observed, expected);
    }

    #[test]
    fn rdk_first_100k_counterexample_matches_rdkit_linear_paths_with_one_bit_per_hash() {
        let observed = observed_active_bits(
            FIRST_100K_COUNTEREXAMPLE,
            RdkFingerprint::default()
                .with_branched_paths(false)
                .with_n_bits_per_hash(1),
        );
        let expected = vec![
            14, 173, 194, 308, 313, 315, 322, 412, 480, 543, 609, 681, 694, 725, 904, 956, 999,
            1057, 1105, 1142, 1147, 1220, 1241, 1308, 1325, 1661, 1780, 2004,
        ];
        assert_eq!(observed, expected);
    }

    #[test]
    fn rdk_first_100k_counterexample_matches_rdkit_without_bond_order() {
        let observed = observed_active_bits(
            FIRST_100K_COUNTEREXAMPLE,
            RdkFingerprint::default()
                .with_n_bits_per_hash(1)
                .with_use_bond_order(false),
        );
        let expected = vec![
            90, 112, 148, 152, 202, 412, 480, 543, 674, 781, 875, 894, 1105, 1233, 1270, 1300,
            1308, 1353, 1582, 1794, 1807, 1809, 1813, 1890, 2019,
        ];
        assert_eq!(observed, expected);
    }

    #[test]
    fn rdk_first_100k_counterexample_matches_rdkit_hypervalent_halogen_bond_codes() {
        let smiles: Smiles = FIRST_100K_COUNTEREXAMPLE
            .parse()
            .expect("counterexample SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);
        let (_, bonds, _) = graph.rdk_fingerprint_graph(true, true);

        let chlorine_atom = (0..graph.atom_count())
            .find(|&atom_id| graph.atom(atom_id).and_then(Atom::element) == Some(Element::Cl))
            .expect("counterexample should contain chlorine");
        let mut chlorine_bond_codes = bonds
            .iter()
            .filter(|bond| bond.source == chlorine_atom || bond.target == chlorine_atom)
            .map(|bond| bond.bond_code)
            .collect::<Vec<_>>();
        chlorine_bond_codes.sort_unstable();

        assert_eq!(chlorine_bond_codes, vec![1, 1, 1, 1]);
    }

    #[test]
    fn rdk_second_500k_counterexample_matches_rdkit_iodine_bond_codes() {
        let smiles: Smiles = SECOND_500K_COUNTEREXAMPLE
            .parse()
            .expect("counterexample SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);
        let (_, bonds, _) = graph.rdk_fingerprint_graph(true, true);

        let iodine_atom = (0..graph.atom_count())
            .find(|&atom_id| graph.atom(atom_id).and_then(Atom::element) == Some(Element::I))
            .expect("counterexample should contain iodine");
        let mut iodine_bond_codes = bonds
            .iter()
            .filter(|bond| bond.source == iodine_atom || bond.target == iodine_atom)
            .map(|bond| bond.bond_code)
            .collect::<Vec<_>>();
        iodine_bond_codes.sort_unstable();

        assert_eq!(iodine_bond_codes, vec![1, 1, 1, 2]);
    }

    #[test]
    fn rdk_third_500k_counterexample_matches_rdkit_iodine_bond_codes() {
        let smiles: Smiles = THIRD_500K_COUNTEREXAMPLE
            .parse()
            .expect("counterexample SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);
        let (_, bonds, _) = graph.rdk_fingerprint_graph(true, true);

        let iodine_atom = (0..graph.atom_count())
            .find(|&atom_id| graph.atom(atom_id).and_then(Atom::element) == Some(Element::I))
            .expect("counterexample should contain iodine");
        let mut iodine_bond_codes = bonds
            .iter()
            .filter(|bond| bond.source == iodine_atom || bond.target == iodine_atom)
            .map(|bond| bond.bond_code)
            .collect::<Vec<_>>();
        iodine_bond_codes.sort_unstable();

        assert_eq!(iodine_bond_codes, vec![1, 1, 1, 2]);
    }

    #[test]
    fn rdk_charge_separated_phosphorylimide_examples_match_rdkit_one_bit_hashes() {
        for (smiles, expected) in [
            ("N=P(=O)O", vec![133, 500, 662, 1030, 1153, 1667, 1799]),
            (
                "CN=P(=O)O",
                vec![469, 500, 628, 669, 1105, 1376, 1667, 1799, 2001],
            ),
            (
                "ON=P(=O)O",
                vec![469, 500, 669, 745, 1039, 1085, 1465, 1667, 1799],
            ),
            ("C=P(=O)C", vec![133, 525, 812, 924, 1803, 1846, 1920]),
            (
                "CC=P(=O)C",
                vec![
                    146, 200, 924, 1208, 1308, 1372, 1667, 1694, 1708, 1803, 1819, 1920,
                ],
            ),
        ] {
            let observed =
                observed_active_bits(smiles, RdkFingerprint::default().with_n_bits_per_hash(1));
            assert_eq!(observed, expected, "failed for {smiles}");
        }
    }

    #[test]
    fn rdk_charge_separated_phosphorylimide_bond_codes_match_rdkit() {
        for (smiles, expected_codes) in [
            ("N=P(=O)O", vec![1, 2, 2]),
            ("CN=P(=O)O", vec![1, 1, 2]),
            ("ON=P(=O)O", vec![1, 1, 2]),
            ("C=P(=O)C", vec![1, 2, 2]),
            ("CC=P(=O)C", vec![1, 1, 2]),
        ] {
            let smiles: Smiles = smiles.parse().expect("fixture SMILES should parse");
            let mut scratch = SmilesRdkitScratch::default();
            let graph = scratch.prepare(&smiles);
            let (_, bonds, _) = graph.rdk_fingerprint_graph(true, true);

            let phosphorus_atom = (0..graph.atom_count())
                .find(|&atom_id| graph.atom(atom_id).and_then(Atom::element) == Some(Element::P))
                .expect("fixture should contain phosphorus");
            let mut phosphorus_bond_codes = bonds
                .iter()
                .filter(|bond| bond.source == phosphorus_atom || bond.target == phosphorus_atom)
                .map(|bond| bond.bond_code)
                .collect::<Vec<_>>();
            phosphorus_bond_codes.sort_unstable();

            assert_eq!(phosphorus_bond_codes, expected_codes, "failed for {smiles}");
        }
    }

    #[test]
    fn rdk_first_1m_counterexample_matches_rdkit_phosphorylimide_bond_codes() {
        let smiles: Smiles = FIRST_1M_COUNTEREXAMPLE
            .parse()
            .expect("counterexample SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);
        let (_, bonds, _) = graph.rdk_fingerprint_graph(true, true);

        let phosphorus_atom = (0..graph.atom_count())
            .find(|&atom_id| graph.atom(atom_id).and_then(Atom::element) == Some(Element::P))
            .expect("counterexample should contain phosphorus");
        let mut phosphorus_bond_codes = bonds
            .iter()
            .filter(|bond| bond.source == phosphorus_atom || bond.target == phosphorus_atom)
            .map(|bond| bond.bond_code)
            .collect::<Vec<_>>();
        phosphorus_bond_codes.sort_unstable();

        assert_eq!(phosphorus_bond_codes, vec![1, 1, 2]);
    }

    #[test]
    fn rdk_cyanophosphoryl_examples_match_rdkit_hashes() {
        for (smiles, expected) in [
            ("C#P=O", vec![133, 1116, 1712]),
            ("NC#P=O", vec![133, 548, 1105, 1116, 1631, 1712]),
            (
                "C(C(=O)O)NC#P=O",
                vec![
                    112, 133, 323, 398, 412, 484, 548, 654, 733, 927, 929, 953, 993, 1039, 1075,
                    1105, 1116, 1154, 1233, 1308, 1339, 1362, 1454, 1537, 1631, 1712, 1721, 1728,
                    1776, 1807, 1901,
                ],
            ),
        ] {
            let observed =
                observed_active_bits(smiles, RdkFingerprint::default().with_n_bits_per_hash(1));
            assert_eq!(observed, expected, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_rdk_matrix_matches_reference_corpus() {
        let fixture = rdkit_rdk_fixture();
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
        assert_eq!(fixture.source.generator, "RDKit RDKFingerprintMol");
        assert_eq!(fixture.source.min_path, 1);
        assert_eq!(fixture.source.max_path, 7);
        assert_eq!(fixture.source.n_bits_per_hash, 2);
        assert!(fixture.source.use_hs);
        assert_eq!(fixture.source.tgt_density, 0.0);
        assert_eq!(fixture.source.min_size, 128);
        assert!(fixture.source.branched_paths);
        assert!(fixture.source.use_bond_order);
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
                let observed_bits = observed_active_bits(smiles, RdkFingerprint::new(case.fp_size));
                assert_eq!(
                    observed_bits, *expected_bits,
                    "failed for fp_size={}, smiles={smiles}",
                    case.fp_size,
                );
            }
        }
    }

    #[test]
    fn rdk_builder_methods_round_trip_custom_settings() {
        let fingerprint = RdkFingerprint::new(513)
            .with_min_path(2)
            .with_max_path(5)
            .with_n_bits_per_hash(3)
            .with_use_hs(false)
            .with_tgt_density(0.42)
            .with_min_size(64)
            .with_branched_paths(false)
            .with_use_bond_order(false);

        assert_eq!(fingerprint.min_path(), 2);
        assert_eq!(fingerprint.max_path(), 5);
        assert_eq!(fingerprint.fp_size(), 513);
        assert_eq!(fingerprint.n_bits_per_hash(), 3);
        assert!(!fingerprint.use_hs());
        assert_eq!(fingerprint.tgt_density(), 0.42);
        assert_eq!(fingerprint.min_size(), 64);
        assert!(!fingerprint.branched_paths());
        assert!(!fingerprint.use_bond_order());
    }

    #[test]
    fn rdk_returns_empty_for_invalid_sizes_or_ranges() {
        let smiles: Smiles = "CCO".parse().expect("fixture SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);

        for fingerprint in [
            RdkFingerprint::new(0),
            RdkFingerprint::new(2048).with_min_path(0),
            RdkFingerprint::new(2048).with_min_path(4).with_max_path(3),
            RdkFingerprint::new(2048).with_n_bits_per_hash(0),
        ] {
            let bits = fingerprint.compute(&graph);
            assert!(bits.active_bits().next().is_none());
            assert_eq!(bits.len(), fingerprint.fp_size());
        }
    }

    #[test]
    fn hashed_rdk_feature_matches_expected_hash_range() {
        let seed = hashed_rdk_feature(
            &[0, 1],
            &[
                RdkFingerprintBond {
                    source: 0,
                    target: 1,
                    bond_code: 1,
                },
                RdkFingerprintBond {
                    source: 1,
                    target: 2,
                    bond_code: 2,
                },
            ],
            &[12, 6, 14],
        );
        assert_eq!(seed, 2_774_995_115);
    }

    #[test]
    fn extra_bits_from_seed_are_deterministic_and_in_range() {
        let mut fingerprint = BitFingerprint::zeros(64);
        set_rdk_feature_bits(&mut fingerprint, 0x1234_5678, 64, 2);
        assert_eq!(fingerprint.active_bits().collect::<Vec<_>>(), vec![3, 56]);
    }
}
