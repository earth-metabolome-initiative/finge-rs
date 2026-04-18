use elements_rs::isotopes::RelativeAtomicMass;
use geometric_traits::traits::{Graph, MonopartiteGraph, MonoplexGraph, TypedNode};
use smiles_parser::{
    atom::{Atom, McesAtomType},
    bond::{Bond, bond_edge::BondEdge},
    smiles::{
        AromaticityAssignmentApplicationError, AromaticityPolicy, RingAtomMembership,
        RingAtomMembershipScratch, Smiles,
    },
};

use crate::traits::{
    AtomPairGraph, EcfpGraph, MolecularAtom, MolecularBond, MolecularGraph, TopologicalTorsionGraph,
};

const ATOM_PAIR_TYPE_BUCKETS: [u8; 15] = [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53];
const ATOM_PAIR_NUM_PI_BITS: u32 = 2;
const ATOM_PAIR_NUM_BRANCH_BITS: u32 = 3;
const ATOM_PAIR_MAX_NUM_PI: u32 = (1 << ATOM_PAIR_NUM_PI_BITS) - 1;
const ATOM_PAIR_MAX_NUM_BRANCHES: u32 = (1 << ATOM_PAIR_NUM_BRANCH_BITS) - 1;

impl MolecularAtom for Atom {
    type AtomType = McesAtomType;

    #[inline]
    fn atom_type(&self) -> Self::AtomType {
        self.node_type()
    }
}

impl MolecularBond for BondEdge {
    type NodeId = usize;
    type BondType = Bond;

    #[inline]
    fn source(&self) -> Self::NodeId {
        self.0
    }

    #[inline]
    fn target(&self) -> Self::NodeId {
        self.1
    }

    #[inline]
    fn bond_type(&self) -> Self::BondType {
        self.2
    }
}

impl MolecularGraph for Smiles {
    type Bond = BondEdge;

    #[inline]
    fn atom(&self, node_id: Self::NodeId) -> Option<&Self::NodeSymbol> {
        self.node_by_id(node_id)
    }

    #[inline]
    fn bonds(&self, node_id: Self::NodeId) -> impl Iterator<Item = Self::Bond> + '_ {
        self.edges_for_node(node_id)
    }
}

#[inline]
fn rdkit_atom_pair_atom_code(smiles: &Smiles, atom_id: usize) -> u32 {
    let Some(atom) = smiles.node_by_id(atom_id) else {
        return 0;
    };

    let branch_count = u32::try_from(smiles.edge_count_for_node(atom_id)).unwrap_or(u32::MAX);
    let num_pi = rdkit_atom_pair_pi_electrons(smiles, atom_id, atom) % ATOM_PAIR_MAX_NUM_PI;
    let type_index = rdkit_atom_pair_type_index(atom);

    (branch_count % ATOM_PAIR_MAX_NUM_BRANCHES)
        | (num_pi << ATOM_PAIR_NUM_BRANCH_BITS)
        | (type_index << (ATOM_PAIR_NUM_BRANCH_BITS + ATOM_PAIR_NUM_PI_BITS))
}

#[inline]
fn rdkit_topological_torsion_atom_code(
    smiles: &Smiles,
    atom_id: usize,
    branch_subtract: u8,
) -> u32 {
    let Some(atom) = smiles.node_by_id(atom_id) else {
        return 0;
    };

    let branch_count = smiles
        .edge_count_for_node(atom_id)
        .saturating_sub(usize::from(branch_subtract));
    let branch_count = u32::try_from(branch_count).unwrap_or(u32::MAX);
    let num_pi = rdkit_atom_pair_pi_electrons(smiles, atom_id, atom) % ATOM_PAIR_MAX_NUM_PI;
    let type_index = rdkit_atom_pair_type_index(atom);

    (branch_count % ATOM_PAIR_MAX_NUM_BRANCHES)
        | (num_pi << ATOM_PAIR_NUM_BRANCH_BITS)
        | (type_index << (ATOM_PAIR_NUM_BRANCH_BITS + ATOM_PAIR_NUM_PI_BITS))
}

#[inline]
fn rdkit_atom_pair_type_index(atom: &Atom) -> u32 {
    let Some(element) = atom.element() else {
        return u32::try_from(ATOM_PAIR_TYPE_BUCKETS.len()).unwrap_or(u32::MAX);
    };

    let atomic_number = u8::from(element);

    ATOM_PAIR_TYPE_BUCKETS
        .iter()
        .position(|&bucket| bucket == atomic_number)
        .map(|index| u32::try_from(index).unwrap_or(u32::MAX))
        .unwrap_or_else(|| u32::try_from(ATOM_PAIR_TYPE_BUCKETS.len()).unwrap_or(u32::MAX))
}

#[inline]
fn hypervalent_halogen_oxoacid_center_double_oxygen_count(
    smiles: &Smiles,
    atom_id: usize,
) -> Option<u32> {
    let atom = smiles.node_by_id(atom_id)?;
    let element = atom.element()?;
    if !matches!(
        element,
        elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I
    ) {
        return None;
    }

    let mut has_single_oxygen = false;
    let mut double_oxygen_count = 0_u32;
    for edge in smiles.edges_for_node(atom_id) {
        let other_id = if edge.source() == atom_id {
            edge.target()
        } else {
            edge.source()
        };
        let Some(other_atom) = smiles.node_by_id(other_id) else {
            continue;
        };
        if other_atom.element() != Some(elements_rs::Element::O) {
            continue;
        }

        match rdkit_bond_order(edge.bond_type()) {
            1 => has_single_oxygen = true,
            2.. => double_oxygen_count += 1,
            _ => {}
        }
    }

    if has_single_oxygen && double_oxygen_count > 0 {
        Some(double_oxygen_count)
    } else {
        None
    }
}

#[inline]
fn normalized_hypervalent_halogen_oxoacid_charge(
    smiles: &Smiles,
    atom_id: usize,
    atom: &Atom,
) -> Option<i8> {
    if let Some(double_oxygen_count) =
        hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, atom_id)
    {
        return i8::try_from(double_oxygen_count).ok();
    }

    if atom.element() == Some(elements_rs::Element::O) {
        for edge in smiles.edges_for_node(atom_id) {
            let other_id = if edge.source() == atom_id {
                edge.target()
            } else {
                edge.source()
            };
            let Some(other_atom) = smiles.node_by_id(other_id) else {
                continue;
            };
            let Some(other_element) = other_atom.element() else {
                continue;
            };
            if !matches!(
                other_element,
                elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I
            ) {
                continue;
            }
            if rdkit_bond_order(edge.bond_type()) >= 2
                && hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, other_id)
                    .is_some()
            {
                return Some(-1);
            }
        }
    }

    None
}

#[inline]
fn rdkit_atom_pair_pi_electrons(smiles: &Smiles, atom_id: usize, atom: &Atom) -> u32 {
    if atom.aromatic() {
        return 1;
    }

    let Some(element) = atom.element() else {
        return 0;
    };

    if element == elements_rs::Element::O {
        for edge in smiles.edges_for_node(atom_id) {
            let other_id = if edge.source() == atom_id {
                edge.target()
            } else {
                edge.source()
            };

            let Some(other_atom) = smiles.node_by_id(other_id) else {
                continue;
            };
            let Some(other_element) = other_atom.element() else {
                continue;
            };
            if !matches!(
                other_element,
                elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I
            ) {
                continue;
            }
            if smiles.edge_count_for_node(other_id) < 2 {
                continue;
            }
            if rdkit_bond_order(edge.bond_type()) >= 2
                && hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, other_id)
                    .is_some()
            {
                return 0;
            }
        }
    }

    let mut degree = 0_usize;
    let mut unsaturation = 0_u32;
    let mut has_single_bond = false;

    for edge in smiles.edges_for_node(atom_id) {
        degree += 1;
        let bond_order = rdkit_bond_order(edge.bond_type());
        unsaturation += bond_order.saturating_sub(1);
        has_single_bond |= bond_order == 1;
    }

    if unsaturation == 0 {
        return 0;
    }
    if !ATOM_PAIR_TYPE_BUCKETS.contains(&u8::from(element)) {
        return rdkit_atom_pair_unknown_element_pi_electrons(
            element,
            degree,
            unsaturation,
            has_single_bond,
        );
    }
    let total_hydrogens =
        usize::from(atom.hydrogen_count()) + usize::from(smiles.implicit_hydrogen_count(atom_id));

    match element {
        elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I => {
            if degree >= 2 { 0 } else { unsaturation }
        }
        elements_rs::Element::P | elements_rs::Element::As | elements_rs::Element::Sb => {
            if (degree <= 2 && !(degree == 2 && total_hydrogens == 2))
                || (degree == 3 && total_hydrogens == 0)
            {
                unsaturation
            } else {
                0
            }
        }
        elements_rs::Element::S | elements_rs::Element::Se | elements_rs::Element::Te => {
            if degree == 2 {
                unsaturation
            } else if unsaturation == 1 {
                if !has_single_bond {
                    1
                } else if atom.charge_value() > 0 {
                    if degree <= 3 { 1 } else { 0 }
                } else if degree == 4 {
                    1
                } else {
                    0
                }
            } else if has_single_bond {
                0
            } else {
                unsaturation
            }
        }
        _ => unsaturation,
    }
}

#[inline]
fn rdkit_atom_invariant(
    smiles: &Smiles,
    atom_id: usize,
    include_ring_membership: bool,
    ring_atom_membership: &RingAtomMembership,
) -> u32 {
    let Some(atom) = smiles.node_by_id(atom_id) else {
        return 0;
    };

    let atomic_number = atom
        .element()
        .map_or(0, |element| u32::from(u8::from(element)));
    let mut explicit_hydrogen_neighbors = 0_u32;
    let mut heavy_degree = 0_u32;
    for edge in smiles.edges_for_node(atom_id) {
        let other_id = if edge.source() == atom_id {
            edge.target()
        } else {
            edge.source()
        };
        let is_hydrogen_neighbor = smiles
            .node_by_id(other_id)
            .and_then(Atom::element)
            .is_some_and(|element| element == elements_rs::Element::H);
        if is_hydrogen_neighbor {
            explicit_hydrogen_neighbors += 1;
        } else {
            heavy_degree += 1;
        }
    }

    let explicit_hydrogens = u32::from(atom.hydrogen_count()) + explicit_hydrogen_neighbors;
    let implicit_hydrogens = u32::from(smiles.implicit_hydrogen_count(atom_id));
    let total_hydrogens = explicit_hydrogens + implicit_hydrogens;
    let total_degree = heavy_degree + total_hydrogens;
    let formal_charge = normalized_hypervalent_halogen_oxoacid_charge(smiles, atom_id, atom)
        .unwrap_or(atom.charge_value()) as i32 as u32;
    let delta_mass = rdkit_delta_mass(atom) as u32;

    if include_ring_membership && ring_atom_membership.contains_atom(atom_id) {
        hash_u32_sequence(&[
            atomic_number,
            total_degree,
            total_hydrogens,
            formal_charge,
            delta_mass,
            1,
        ])
    } else {
        hash_u32_sequence(&[
            atomic_number,
            total_degree,
            total_hydrogens,
            formal_charge,
            delta_mass,
        ])
    }
}

#[inline]
fn rdkit_delta_mass(atom: &Atom) -> i32 {
    let Some(element) = atom.element() else {
        return 0;
    };

    let standard_atomic_mass = rdkit_atomic_weight(element);
    let atom_mass = atom
        .isotope_mass_number()
        .and_then(|_| atom.isotope().ok())
        .map(|isotope| isotope.relative_atomic_mass())
        .unwrap_or(standard_atomic_mass);

    (atom_mass - standard_atomic_mass) as i32
}

#[inline]
fn rdkit_atomic_weight(element: elements_rs::Element) -> f64 {
    match element {
        elements_rs::Element::Tc => 98.0,
        elements_rs::Element::Pm => 145.0,
        elements_rs::Element::Po => 209.0,
        elements_rs::Element::At => 210.0,
        elements_rs::Element::Rn => 222.0,
        elements_rs::Element::Fr => 223.0,
        elements_rs::Element::Ra => 226.0,
        elements_rs::Element::Ac => 227.0,
        elements_rs::Element::Np => 237.0,
        elements_rs::Element::Pu => 244.0,
        elements_rs::Element::Am => 243.0,
        elements_rs::Element::Cm => 247.0,
        elements_rs::Element::Bk => 247.0,
        elements_rs::Element::Cf => 251.0,
        elements_rs::Element::Es => 252.0,
        elements_rs::Element::Fm => 257.0,
        elements_rs::Element::Md => 258.0,
        elements_rs::Element::No => 259.0,
        elements_rs::Element::Lr => 262.0,
        _ => (element.standard_atomic_weight() * 1000.0).round() / 1000.0,
    }
}

#[inline]
fn hash_u32_sequence(values: &[u32]) -> u32 {
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

/// Error raised while preparing a `smiles-parser` graph for RDKit-parity
/// fingerprinting.
#[derive(Debug)]
pub enum SmilesPreparationError {
    /// Applying the RDKit-default aromaticity assignment failed.
    Aromaticity(AromaticityAssignmentApplicationError),
}

impl core::fmt::Display for SmilesPreparationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Aromaticity(error) => {
                write!(
                    f,
                    "failed to apply RDKit-default aromaticity for fingerprinting: {error}"
                )
            }
        }
    }
}

fn normalize_smiles_for_rdkit(smiles: &Smiles) -> Result<Smiles, SmilesPreparationError> {
    let has_aromatic_labels = smiles.nodes().iter().any(Atom::aromatic)
        || (0..smiles.atom_count()).any(|atom_id| {
            smiles
                .edges_for_node(atom_id)
                .any(|bond| bond.bond_type() == Bond::Aromatic)
        });

    let base = if has_aromatic_labels {
        if let Ok(kekule) = smiles.kekulize_standalone() {
            kekule
        } else {
            return Ok(smiles.clone());
        }
    } else {
        smiles.clone()
    };

    let perception = base
        .perceive_aromaticity_for(AromaticityPolicy::RdkitDefault)
        .map_err(SmilesPreparationError::Aromaticity)?;
    Ok(perception.into_aromaticized())
}

/// Reusable scratch for batch RDKit-parity processing over `smiles-parser`
/// graphs.
///
/// This keeps the atom-only ring-membership output and DFS scratch buffers
/// alive across molecules so callers can avoid reallocating them on every
/// preparation step.
#[derive(Debug, Default)]
pub struct SmilesRdkitScratch {
    normalized_smiles: Option<Smiles>,
    ring_atom_membership: RingAtomMembership,
    ring_atom_membership_scratch: RingAtomMembershipScratch,
}

impl SmilesRdkitScratch {
    /// Prepares a `smiles-parser` graph for RDKit-parity fingerprinting while
    /// reusing the normalized graph and atom-ring-membership buffers held by
    /// this scratch object.
    #[inline]
    #[must_use]
    pub fn prepare<'a>(&'a mut self, smiles: &Smiles) -> SmilesRdkitGraph<'a> {
        self.try_prepare(smiles)
            .expect("fingerprint preparation should succeed for valid SMILES inputs")
    }

    /// Prepares a `smiles-parser` graph for RDKit-parity fingerprinting.
    ///
    /// This first normalizes aromaticity to the RDKit-default model, then
    /// populates reusable ring-membership buffers for ECFP.
    pub fn try_prepare<'a>(
        &'a mut self,
        smiles: &Smiles,
    ) -> Result<SmilesRdkitGraph<'a>, SmilesPreparationError> {
        self.normalized_smiles = Some(normalize_smiles_for_rdkit(smiles)?);
        let normalized = self
            .normalized_smiles
            .as_ref()
            .expect("normalized graph should be present after successful preparation");

        normalized.write_ring_atom_membership(
            &mut self.ring_atom_membership,
            &mut self.ring_atom_membership_scratch,
        );

        Ok(SmilesRdkitGraph {
            smiles: normalized,
            ring_atom_membership: &self.ring_atom_membership,
        })
    }
}

/// Scratch-prepared RDKit-normalized `smiles-parser` graph view.
///
/// Create this through [`SmilesRdkitScratch::prepare`] and pass it to the
/// normal [`Fingerprint::compute`](crate::Fingerprint::compute) path.
#[derive(Debug, Clone, Copy)]
pub struct SmilesRdkitGraph<'a> {
    smiles: &'a Smiles,
    ring_atom_membership: &'a RingAtomMembership,
}

impl Graph for SmilesRdkitGraph<'_> {
    #[inline]
    fn has_nodes(&self) -> bool {
        self.smiles.has_nodes()
    }

    #[inline]
    fn has_edges(&self) -> bool {
        self.smiles.has_edges()
    }
}

impl MonoplexGraph for SmilesRdkitGraph<'_> {
    type Edge = <Smiles as MonoplexGraph>::Edge;
    type Edges = <Smiles as MonoplexGraph>::Edges;

    #[inline]
    fn edges(&self) -> &Self::Edges {
        self.smiles.edges()
    }
}

impl MonopartiteGraph for SmilesRdkitGraph<'_> {
    type NodeId = <Smiles as MonopartiteGraph>::NodeId;
    type NodeSymbol = <Smiles as MonopartiteGraph>::NodeSymbol;
    type Nodes = <Smiles as MonopartiteGraph>::Nodes;

    #[inline]
    fn nodes_vocabulary(&self) -> &Self::Nodes {
        self.smiles.nodes_vocabulary()
    }
}

impl MolecularGraph for SmilesRdkitGraph<'_> {
    type Bond = BondEdge;

    #[inline]
    fn atom(&self, node_id: Self::NodeId) -> Option<&Self::NodeSymbol> {
        self.smiles.node_by_id(node_id)
    }

    #[inline]
    fn bonds(&self, node_id: Self::NodeId) -> impl Iterator<Item = Self::Bond> + '_ {
        self.smiles.edges_for_node(node_id)
    }
}

impl EcfpGraph for SmilesRdkitGraph<'_> {
    #[inline]
    fn ecfp_atom_invariant(&self, atom_id: usize, include_ring_membership: bool) -> u32 {
        rdkit_atom_invariant(
            self.smiles,
            atom_id,
            include_ring_membership,
            self.ring_atom_membership,
        )
    }

    #[inline]
    fn ecfp_bond_invariant(&self, bond: &Self::Bond, use_bond_types: bool) -> u32 {
        if !use_bond_types {
            return 1;
        }

        let begin_atom_id = bond.source();
        let end_atom_id = bond.target();
        let begin_atom = self.smiles.node_by_id(begin_atom_id);
        let end_atom = self.smiles.node_by_id(end_atom_id);
        let is_halogen_oxygen_bond = matches!(
            (
                begin_atom.and_then(Atom::element),
                end_atom.and_then(Atom::element)
            ),
            (
                Some(elements_rs::Element::O),
                Some(elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I),
            ) | (
                Some(elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I),
                Some(elements_rs::Element::O),
            )
        );
        let is_hypervalent_oxygen_halogen_bond = is_halogen_oxygen_bond
            && (hypervalent_halogen_oxoacid_center_double_oxygen_count(self.smiles, begin_atom_id)
                .is_some()
                || hypervalent_halogen_oxoacid_center_double_oxygen_count(
                    self.smiles,
                    end_atom_id,
                )
                .is_some());
        if is_hypervalent_oxygen_halogen_bond {
            return 1;
        }

        match bond.bond_type() {
            Bond::Single | Bond::Up | Bond::Down => 1,
            Bond::Double => 2,
            Bond::Triple => 3,
            Bond::Quadruple => 4,
            Bond::Aromatic => 12,
        }
    }
}

impl AtomPairGraph for Smiles {
    #[inline]
    fn atom_pair_atom_code(&self, atom_id: usize) -> u32 {
        rdkit_atom_pair_atom_code(self, atom_id)
    }
}

impl AtomPairGraph for SmilesRdkitGraph<'_> {
    #[inline]
    fn atom_pair_atom_code(&self, atom_id: usize) -> u32 {
        rdkit_atom_pair_atom_code(self.smiles, atom_id)
    }
}

impl TopologicalTorsionGraph for Smiles {
    #[inline]
    fn topological_torsion_atom_code(&self, atom_id: usize, branch_subtract: u8) -> u32 {
        rdkit_topological_torsion_atom_code(self, atom_id, branch_subtract)
    }
}

impl TopologicalTorsionGraph for SmilesRdkitGraph<'_> {
    #[inline]
    fn topological_torsion_atom_code(&self, atom_id: usize, branch_subtract: u8) -> u32 {
        rdkit_topological_torsion_atom_code(self.smiles, atom_id, branch_subtract)
    }
}

#[inline]
fn rdkit_bond_order(bond: Bond) -> u32 {
    match bond {
        Bond::Single | Bond::Up | Bond::Down | Bond::Aromatic => 1,
        Bond::Double => 2,
        Bond::Triple => 3,
        Bond::Quadruple => 4,
    }
}

#[inline]
fn rdkit_atom_pair_unknown_element_pi_electrons(
    element: elements_rs::Element,
    degree: usize,
    unsaturation: u32,
    has_single_bond: bool,
) -> u32 {
    match element {
        // Group 7: terminal multiple-bond cases contribute their unsaturation,
        // degree-2 all-multiple-bond cases carry unsaturation except for
        // dicarbyne-like patterns, which RDKit buckets as one pi electron,
        // and mixed single/multiple coordination collapses to zero.
        elements_rs::Element::Mn
        | elements_rs::Element::Tc
        | elements_rs::Element::Re => match degree {
            1 => unsaturation,
            2 if !has_single_bond => {
                if unsaturation == 4 {
                    1
                } else {
                    unsaturation
                }
            }
            _ => 0,
        },
        // Group 6: terminal and degree-2 oxo-like cases contribute, but degree-3
        // mono-oxo hydroxy patterns and degree-4 dioxo dihydroxy patterns do not.
        elements_rs::Element::Cr
        | elements_rs::Element::Mo
        | elements_rs::Element::W
        // Early actinyl-forming actinides follow the same RDKit wildcard rule.
        | elements_rs::Element::Th
        | elements_rs::Element::U
        | elements_rs::Element::Np => match degree {
            0 => 0,
            1 | 2 => unsaturation,
            3 if unsaturation >= 2 => unsaturation,
            3 => 0,
            _ => 0,
        },
        // Group 5: unsaturation is carried through up to degree 3, but
        // degree-4 dioxo dihydroxy patterns collapse to zero.
        elements_rs::Element::V
        | elements_rs::Element::Nb
        | elements_rs::Element::Ta => {
            if degree <= 3 {
                unsaturation
            } else {
                0
            }
        }
        // Group 8: terminal and degree-2 oxo-like cases collapse to zero, the
        // degree-3 mono-oxo hydroxy pattern contributes one, and degree-4 dioxo
        // dihydroxy patterns contribute two.
        elements_rs::Element::Fe
        | elements_rs::Element::Ru
        | elements_rs::Element::Os => match degree {
            0 => 0,
            1 if unsaturation >= 2 => unsaturation,
            1 => 0,
            2 => 0,
            3 if unsaturation == 1 => 1,
            3 => 0,
            4 if unsaturation == 2 => 2,
            4 => 0,
            _ => 0,
        },
        // Group 9: terminal oxo-like cases collapse to zero, but higher-degree
        // mixed oxo/hydroxy patterns contribute one or two pi electrons in the
        // same way RDKit buckets cobalt, rhodium, and iridium.
        elements_rs::Element::Co
        | elements_rs::Element::Rh
        | elements_rs::Element::Ir => match degree {
            0 | 1 => 0,
            2 if unsaturation == 1 => 1,
            2 => 0,
            3 if unsaturation == 1 => 1,
            3 => 0,
            4 if unsaturation == 2 => 2,
            4 => 0,
            _ => 0,
        },
        // Group 10: only all-double-bond patterns contribute; terminal and
        // multi-terminal triple-bond patterns collapse to zero.
        elements_rs::Element::Ni
        | elements_rs::Element::Pd
        | elements_rs::Element::Pt
        // Dysprosium follows the same RDKit wildcard rule for the cases seen
        // in the parity corpus.
        | elements_rs::Element::Dy => {
            if unsaturation == degree as u32 {
                unsaturation
            } else {
                0
            }
        }
        // Terbium collapses to zero in the RDKit wildcard cases observed so
        // far, including terminal multiple bonds.
        elements_rs::Element::Tb => 0,
        // Other wildcard-type elements follow the unsaturation count directly.
        _ => unsaturation,
    }
}

#[cfg(test)]
mod tests {
    use alloc::{vec, vec::Vec};

    use smiles_parser::smiles::Smiles;

    use super::SmilesRdkitScratch;
    use crate::{EcfpGraph, MolecularGraph};

    fn observed_connectivity_invariants(smiles: &str, include_ring_membership: bool) -> Vec<u32> {
        let smiles: Smiles = smiles.parse().expect("fixture SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let graph = scratch.prepare(&smiles);

        (0..graph.atom_count())
            .map(|atom_id| graph.ecfp_atom_invariant(atom_id, include_ring_membership))
            .collect()
    }

    #[test]
    fn rdkit_connectivity_invariants_match_reference_fixtures() {
        for (smiles, expected) in [
            ("CN", vec![2246728737, 847957139]),
            ("O=C=O", vec![864942730, 2245900962, 864942730]),
            (
                "c1ccncc1",
                vec![
                    3218693969, 3218693969, 3218693969, 2041434490, 3218693969, 3218693969,
                ],
            ),
            ("[NH4+]", vec![847680145]),
            ("[13CH4]", vec![2246733040]),
            ("C[NH3+]", vec![2246728737, 847694221]),
        ] {
            let observed = observed_connectivity_invariants(smiles, true);
            assert_eq!(observed, expected, "failed for {smiles}");
        }
    }

    #[test]
    fn rdkit_connectivity_invariants_without_ring_membership_match_reference_fixtures() {
        for (smiles, expected) in [
            (
                "c1ccccc1",
                vec![
                    2246703798, 2246703798, 2246703798, 2246703798, 2246703798, 2246703798,
                ],
            ),
            (
                "C1CCCCC1",
                vec![
                    2245384272, 2245384272, 2245384272, 2245384272, 2245384272, 2245384272,
                ],
            ),
            (
                "C1=CC=CN=C1",
                vec![
                    2246703798, 2246703798, 2246703798, 2246703798, 847336149, 2246703798,
                ],
            ),
            (
                "OC1CCCCC1",
                vec![
                    864662311, 2245273601, 2245384272, 2245384272, 2245384272, 2245384272,
                    2245384272,
                ],
            ),
        ] {
            let observed = observed_connectivity_invariants(smiles, false);
            assert_eq!(observed, expected, "failed for {smiles}");
        }
    }
}
