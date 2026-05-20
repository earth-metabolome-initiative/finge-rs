use alloc::{vec, vec::Vec};

use elements_rs::{AllowedValences, ValenceElectrons, isotopes::RelativeAtomicMass};
use geometric_traits::traits::{Graph, MonopartiteGraph, MonoplexGraph, TypedNode};
use smiles_parser::{
    atom::{Atom, McesAtomType},
    bond::{Bond, bond_edge::BondEdge},
    smiles::{
        AromaticityAssignmentApplicationError, AromaticityPolicy, RingAtomMembership,
        RingAtomMembershipScratch, Smiles,
    },
};

use crate::{
    atom_invariant_fields::AtomInvariantFields,
    traits::{
        AtomPairGraph, EcfpGraph, MolecularAtom, MolecularBond, MolecularGraph,
        RdkFingerprintGraph, TopologicalTorsionGraph,
    },
};

const ATOM_PAIR_TYPE_BUCKETS: [u8; 15] = [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53];
const ATOM_PAIR_NUM_PI_BITS: u32 = 2;
const ATOM_PAIR_NUM_BRANCH_BITS: u32 = 3;
const ATOM_PAIR_NUM_TYPE_BITS: u32 = 4;
const ATOM_PAIR_CODE_SIZE: u32 =
    ATOM_PAIR_NUM_BRANCH_BITS + ATOM_PAIR_NUM_PI_BITS + ATOM_PAIR_NUM_TYPE_BITS;
const ATOM_PAIR_MAX_NUM_PI: u32 = (1 << ATOM_PAIR_NUM_PI_BITS) - 1;
const ATOM_PAIR_MAX_NUM_BRANCHES: u32 = (1 << ATOM_PAIR_NUM_BRANCH_BITS) - 1;
const TOPOLOGICAL_TORSION_CODE_MODULUS: u32 = (1 << ATOM_PAIR_CODE_SIZE) - 1;

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
        BondEdge::source(*self)
    }

    #[inline]
    fn target(&self) -> Self::NodeId {
        BondEdge::target(*self)
    }

    #[inline]
    fn bond_type(&self) -> Self::BondType {
        BondEdge::bond_type(*self)
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
    rdkit_atom_pair_atom_code_with_options(smiles, atom_id, false)
}

#[inline]
fn rdkit_atom_pair_atom_code_with_options(
    smiles: &Smiles,
    atom_id: usize,
    collapse_plain_explicit_hydrogens: bool,
) -> u32 {
    let Some(atom) = smiles.node_by_id(atom_id) else {
        return 0;
    };

    let branch_count =
        rdkit_atom_pair_branch_count(smiles, atom_id, collapse_plain_explicit_hydrogens);
    let num_pi =
        rdkit_atom_pair_pi_electrons(smiles, atom_id, atom, collapse_plain_explicit_hydrogens)
            % ATOM_PAIR_MAX_NUM_PI;
    let type_index = rdkit_atom_pair_type_index(atom);

    (branch_count % ATOM_PAIR_MAX_NUM_BRANCHES)
        | (num_pi << ATOM_PAIR_NUM_BRANCH_BITS)
        | (type_index << (ATOM_PAIR_NUM_BRANCH_BITS + ATOM_PAIR_NUM_PI_BITS))
}

#[inline]
fn rdkit_atom_pair_branch_count(
    smiles: &Smiles,
    atom_id: usize,
    collapse_plain_explicit_hydrogens: bool,
) -> u32 {
    let branch_count = smiles
        .edges_for_node(atom_id)
        .filter(|edge| {
            let other_id = if edge.source() == atom_id {
                edge.target()
            } else {
                edge.source()
            };
            !collapse_plain_explicit_hydrogens
                || !rdkit_rdk_atom_is_collapsible_hydrogen(smiles, other_id)
        })
        .count();

    u32::try_from(branch_count).unwrap_or(u32::MAX)
}

#[inline]
fn rdkit_topological_torsion_atom_code(
    smiles: &Smiles,
    atom_id: usize,
    branch_subtract: u8,
) -> u32 {
    rdkit_topological_torsion_atom_code_with_options(smiles, atom_id, branch_subtract, false)
}

#[inline]
fn rdkit_topological_torsion_atom_code_with_options(
    smiles: &Smiles,
    atom_id: usize,
    branch_subtract: u8,
    collapse_plain_explicit_hydrogens: bool,
) -> u32 {
    if smiles.node_by_id(atom_id).is_none() {
        return 0;
    }

    let atom_code =
        rdkit_atom_pair_atom_code_with_options(smiles, atom_id, collapse_plain_explicit_hydrogens);
    // RDKit computes TT path codes from `getAtomCode(atom, 0) - 2` as an
    // unsigned invariant before folding modulo `2^codeSize - 1`.
    let code = atom_code.wrapping_sub(2) % TOPOLOGICAL_TORSION_CODE_MODULUS + 1;
    code.saturating_sub(u32::from(branch_subtract.saturating_sub(1)))
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

    let mut single_oxygen_count = 0_u32;
    let mut double_oxygen_count = 0_u32;
    let mut non_oxygen_ligand_count = 0_u32;
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
            non_oxygen_ligand_count += 1;
            continue;
        }

        match rdkit_bond_order(edge.bond_type()) {
            1 => single_oxygen_count += 1,
            2.. => double_oxygen_count += 1,
            _ => {}
        }
    }

    match element {
        elements_rs::Element::Cl | elements_rs::Element::Br => {
            (single_oxygen_count > 0 && double_oxygen_count > 0).then_some(double_oxygen_count)
        }
        elements_rs::Element::I => {
            let keep_single_double = double_oxygen_count == 1
                && (single_oxygen_count == 2 || non_oxygen_ligand_count > 0);
            (single_oxygen_count > 0 && double_oxygen_count > 0 && !keep_single_double)
                .then_some(double_oxygen_count)
        }
        _ => None,
    }
}

#[inline]
fn hydrogenated_hypervalent_halogen_center_double_oxygen_count(
    smiles: &Smiles,
    atom_id: usize,
) -> Option<u32> {
    let atom = smiles.node_by_id(atom_id)?;
    let element = atom.element()?;
    if !matches!(
        element,
        elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I
    ) || atom.hydrogen_count() == 0
    {
        return None;
    }

    let mut single_oxygen_count = 0_u32;
    let mut double_oxygen_count = 0_u32;
    let mut non_oxygen_ligand_count = 0_u32;
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
            non_oxygen_ligand_count += 1;
            continue;
        }

        match rdkit_bond_order(edge.bond_type()) {
            1 => single_oxygen_count += 1,
            2.. => double_oxygen_count += 1,
            _ => {}
        }
    }

    let charge_separated = match element {
        elements_rs::Element::Cl | elements_rs::Element::Br => {
            double_oxygen_count > 0 && non_oxygen_ligand_count == 0
        }
        elements_rs::Element::I => {
            double_oxygen_count > 0 && single_oxygen_count >= 2 && non_oxygen_ligand_count == 0
        }
        _ => false,
    };

    charge_separated.then_some(double_oxygen_count)
}

#[inline]
fn rdkit_atom_pair_valence_bond_order(smiles: &Smiles, bond: BondEdge) -> u32 {
    if is_hypervalent_halogen_oxygen_bond(smiles, bond)
        || is_charge_separated_phosphorylimide_oxygen_bond(smiles, bond)
    {
        return 1;
    }

    rdkit_bond_order(bond.bond_type())
}

#[inline]
fn charge_separated_phosphoryl_center_double_oxygen_count(
    smiles: &Smiles,
    atom_id: usize,
) -> Option<u32> {
    let atom = smiles.node_by_id(atom_id)?;
    if atom.element() != Some(elements_rs::Element::P) {
        return None;
    }
    if atom.charge_value() != 0 {
        return None;
    }

    let mut double_oxygen_count = 0_u32;
    let mut has_charge_separating_pi_ligand = false;
    for edge in smiles.edges_for_node(atom_id) {
        let other_id = if edge.source() == atom_id {
            edge.target()
        } else {
            edge.source()
        };
        let Some(other_atom) = smiles.node_by_id(other_id) else {
            continue;
        };
        let bond_order = rdkit_bond_order(edge.bond_type());
        match other_atom.element() {
            Some(elements_rs::Element::O) if bond_order >= 2 && !edge.is_aromatic() => {
                double_oxygen_count += 1;
            }
            Some(elements_rs::Element::C | elements_rs::Element::N)
                if bond_order == 2 && !edge.is_aromatic() =>
            {
                let has_non_hydrogen_substituent =
                    smiles.edges_for_node(other_id).any(|other_edge| {
                        let neighbor_id = if other_edge.source() == other_id {
                            other_edge.target()
                        } else {
                            other_edge.source()
                        };
                        neighbor_id != atom_id
                    });
                if has_non_hydrogen_substituent {
                    has_charge_separating_pi_ligand = true;
                }
            }
            _ => {}
        }
    }

    (double_oxygen_count > 0 && (atom.aromatic() || has_charge_separating_pi_ligand))
        .then_some(double_oxygen_count)
}

#[inline]
fn normalized_rdkit_charge(smiles: &Smiles, atom_id: usize, atom: &Atom) -> Option<i8> {
    if let Some(double_oxygen_count) =
        hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, atom_id)
    {
        return i8::try_from(double_oxygen_count).ok();
    }
    if let Some(double_oxygen_count) =
        hydrogenated_hypervalent_halogen_center_double_oxygen_count(smiles, atom_id)
    {
        return i8::try_from(double_oxygen_count).ok();
    }

    if let Some(double_oxygen_count) =
        charge_separated_phosphoryl_center_double_oxygen_count(smiles, atom_id)
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
            if matches!(
                other_element,
                elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I
            ) {
                if rdkit_bond_order(edge.bond_type()) >= 2
                    && hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, other_id)
                        .is_some()
                {
                    return Some(-1);
                }
                if rdkit_bond_order(edge.bond_type()) >= 2
                    && hydrogenated_hypervalent_halogen_center_double_oxygen_count(smiles, other_id)
                        .is_some()
                {
                    return Some(-1);
                }
            }
            if other_element == elements_rs::Element::P
                && rdkit_bond_order(edge.bond_type()) >= 2
                && !edge.is_aromatic()
                && charge_separated_phosphoryl_center_double_oxygen_count(smiles, other_id)
                    .is_some()
            {
                return Some(-1);
            }
        }
    }

    None
}

#[derive(Debug, Clone, Copy)]
struct RdkitAtomPairValenceState {
    explicit_valence: u32,
    physical_bonds: u32,
    graph_degree: u32,
    total_degree: u32,
    total_valence: u32,
    formal_charge: i8,
}

#[inline]
fn rdkit_atom_pair_valence_state(
    smiles: &Smiles,
    atom_id: usize,
    atom: &Atom,
    collapse_plain_explicit_hydrogens: bool,
) -> RdkitAtomPairValenceState {
    let formal_charge =
        normalized_rdkit_charge(smiles, atom_id, atom).unwrap_or(atom.charge_value());
    let collapsed_hydrogens = if collapse_plain_explicit_hydrogens {
        smiles
            .edges_for_node(atom_id)
            .filter(|edge| {
                let other_id = if edge.source() == atom_id {
                    edge.target()
                } else {
                    edge.source()
                };
                rdkit_rdk_atom_is_collapsible_hydrogen(smiles, other_id)
            })
            .count()
    } else {
        0
    };
    let explicit_hydrogens =
        u32::from(atom.hydrogen_count()) + u32::try_from(collapsed_hydrogens).unwrap_or(u32::MAX);
    let implicit_hydrogens = u32::from(smiles.implicit_hydrogen_count(atom_id));
    let mut explicit_valence = explicit_hydrogens;
    let mut physical_bonds = explicit_hydrogens;
    let mut graph_degree = 0_u32;

    for edge in smiles.edges_for_node(atom_id) {
        let other_id = if edge.source() == atom_id {
            edge.target()
        } else {
            edge.source()
        };
        if collapse_plain_explicit_hydrogens
            && rdkit_rdk_atom_is_collapsible_hydrogen(smiles, other_id)
        {
            continue;
        }

        let bond_order = rdkit_atom_pair_valence_bond_order(smiles, edge);
        explicit_valence = explicit_valence.saturating_add(bond_order);
        if bond_order != 0 {
            physical_bonds = physical_bonds.saturating_add(1);
        }
        graph_degree = graph_degree.saturating_add(1);
    }

    RdkitAtomPairValenceState {
        explicit_valence,
        physical_bonds,
        graph_degree,
        total_degree: graph_degree
            .saturating_add(explicit_hydrogens)
            .saturating_add(implicit_hydrogens),
        total_valence: explicit_valence.saturating_add(implicit_hydrogens),
        formal_charge,
    }
}

#[inline]
fn rdkit_atom_pair_outer_electrons(element: elements_rs::Element) -> i32 {
    // RDKit's AtomPair/TopologicalTorsion atom code depends on the
    // PeriodicTable::getNOuterElecs() compatibility table. That table is not
    // the same contract as elements-rs::ValenceElectrons, which is a
    // group-style chemistry descriptor.
    let rdkit_outer_electrons = match element {
        elements_rs::Element::Zn | elements_rs::Element::Cd | elements_rs::Element::Hg => 2,
        elements_rs::Element::Pr | elements_rs::Element::Pa => 3,
        elements_rs::Element::Nd | elements_rs::Element::U => 4,
        elements_rs::Element::Pm | elements_rs::Element::Np => 5,
        elements_rs::Element::Sm | elements_rs::Element::Pu => 6,
        elements_rs::Element::Eu | elements_rs::Element::Am => 7,
        elements_rs::Element::Gd | elements_rs::Element::Cm => 8,
        elements_rs::Element::Tb | elements_rs::Element::Bk => 9,
        elements_rs::Element::Dy | elements_rs::Element::Cf => 10,
        elements_rs::Element::Ho | elements_rs::Element::Es => 11,
        elements_rs::Element::Er | elements_rs::Element::Fm => 12,
        elements_rs::Element::Tm | elements_rs::Element::Md => 13,
        elements_rs::Element::Yb | elements_rs::Element::No => 14,
        elements_rs::Element::Lu | elements_rs::Element::Lr => 15,
        elements_rs::Element::Rf
        | elements_rs::Element::Db
        | elements_rs::Element::Sg
        | elements_rs::Element::Bh
        | elements_rs::Element::Hs
        | elements_rs::Element::Mt
        | elements_rs::Element::Ds
        | elements_rs::Element::Rg
        | elements_rs::Element::Cn
        | elements_rs::Element::Nh
        | elements_rs::Element::Fl
        | elements_rs::Element::Mc
        | elements_rs::Element::Lv
        | elements_rs::Element::Ts
        | elements_rs::Element::Og => 2,
        _ => element.valence_electrons(),
    };

    i32::from(rdkit_outer_electrons)
}

#[inline]
fn rdkit_atom_pair_default_valence(element: elements_rs::Element) -> i32 {
    // RDKit uses -1 as an open/unspecified valence sentinel for many metals.
    // Keep that behavior local to the RDKit-compatible atom-code path.
    if matches!(
        element,
        elements_rs::Element::Sc
            | elements_rs::Element::Ti
            | elements_rs::Element::V
            | elements_rs::Element::Cr
            | elements_rs::Element::Mn
            | elements_rs::Element::Fe
            | elements_rs::Element::Co
            | elements_rs::Element::Ni
            | elements_rs::Element::Cu
            | elements_rs::Element::Zn
            | elements_rs::Element::Y
            | elements_rs::Element::Zr
            | elements_rs::Element::Nb
            | elements_rs::Element::Mo
            | elements_rs::Element::Tc
            | elements_rs::Element::Ru
            | elements_rs::Element::Rh
            | elements_rs::Element::Pd
            | elements_rs::Element::Ag
            | elements_rs::Element::Cd
            | elements_rs::Element::La
            | elements_rs::Element::Ce
            | elements_rs::Element::Pr
            | elements_rs::Element::Nd
            | elements_rs::Element::Pm
            | elements_rs::Element::Sm
            | elements_rs::Element::Eu
            | elements_rs::Element::Gd
            | elements_rs::Element::Tb
            | elements_rs::Element::Dy
            | elements_rs::Element::Ho
            | elements_rs::Element::Er
            | elements_rs::Element::Tm
            | elements_rs::Element::Yb
            | elements_rs::Element::Lu
            | elements_rs::Element::Hf
            | elements_rs::Element::Ta
            | elements_rs::Element::W
            | elements_rs::Element::Re
            | elements_rs::Element::Os
            | elements_rs::Element::Ir
            | elements_rs::Element::Pt
            | elements_rs::Element::Au
            | elements_rs::Element::Hg
            | elements_rs::Element::Tl
            | elements_rs::Element::Ac
            | elements_rs::Element::Th
            | elements_rs::Element::Pa
            | elements_rs::Element::U
            | elements_rs::Element::Np
            | elements_rs::Element::Pu
            | elements_rs::Element::Am
            | elements_rs::Element::Cm
            | elements_rs::Element::Bk
            | elements_rs::Element::Cf
            | elements_rs::Element::Es
            | elements_rs::Element::Fm
            | elements_rs::Element::Md
            | elements_rs::Element::No
    ) {
        return -1;
    }

    element
        .allowed_valences()
        .first()
        .copied()
        .map_or(-1, i32::from)
}

#[inline]
fn rdkit_atom_pair_num_bonds_plus_lone_pairs(
    element: elements_rs::Element,
    state: RdkitAtomPairValenceState,
) -> i32 {
    let degree = i32::try_from(state.total_degree).unwrap_or(i32::MAX);
    if u8::from(element) <= 1 {
        return degree;
    }

    let outer_electrons = rdkit_atom_pair_outer_electrons(element);
    let total_valence = i32::try_from(state.total_valence).unwrap_or(i32::MAX);
    let formal_charge = i32::from(state.formal_charge);
    let free_electrons = outer_electrons.saturating_sub(total_valence + formal_charge);

    degree + free_electrons / 2
}

#[inline]
fn rdkit_atom_pair_count_atom_electrons(
    element: elements_rs::Element,
    state: RdkitAtomPairValenceState,
) -> i32 {
    let default_valence = rdkit_atom_pair_default_valence(element);
    if default_valence <= 1 {
        return -1;
    }

    let degree = i32::try_from(state.total_degree).unwrap_or(i32::MAX);
    if degree > 3 {
        return -1;
    }

    let lone_pair_electrons = (rdkit_atom_pair_outer_electrons(element) - default_valence)
        .saturating_sub(i32::from(state.formal_charge))
        .max(0);
    let mut electrons = (default_valence - degree) + lone_pair_electrons;

    if electrons > 1 {
        let explicit_valence = i32::try_from(state.explicit_valence).unwrap_or(i32::MAX);
        let graph_degree = i32::try_from(state.graph_degree).unwrap_or(i32::MAX);
        if explicit_valence.saturating_sub(graph_degree) > 1 {
            electrons = 1;
        }
    }

    electrons
}

#[inline]
fn rdkit_atom_pair_is_conjugation_candidate(
    element: elements_rs::Element,
    state: RdkitAtomPairValenceState,
) -> bool {
    let default_valence = rdkit_atom_pair_default_valence(element);
    if state.formal_charge == 0
        && default_valence >= 0
        && state.total_valence > u32::try_from(default_valence).unwrap_or(u32::MAX)
    {
        return false;
    }

    let outer_electrons = rdkit_atom_pair_outer_electrons(element);
    (u8::from(element) <= 10
        || (outer_electrons != 5 && outer_electrons != 6)
        || (outer_electrons == 6 && state.total_degree < 2))
        && rdkit_atom_pair_count_atom_electrons(element, state) > 0
}

#[inline]
fn rdkit_atom_pair_neighbor_id(edge: BondEdge, atom_id: usize) -> usize {
    if edge.source() == atom_id {
        edge.target()
    } else {
        edge.source()
    }
}

#[inline]
fn rdkit_atom_pair_conjugation_candidate_state(
    smiles: &Smiles,
    atom_id: usize,
    collapse_plain_explicit_hydrogens: bool,
) -> Option<RdkitAtomPairValenceState> {
    let atom = smiles.node_by_id(atom_id)?;
    let element = atom.element()?;
    let state =
        rdkit_atom_pair_valence_state(smiles, atom_id, atom, collapse_plain_explicit_hydrogens);
    rdkit_atom_pair_is_conjugation_candidate(element, state).then_some(state)
}

#[inline]
fn rdkit_atom_pair_conjugated_by_center(
    smiles: &Smiles,
    center_id: usize,
    marked_neighbor_id: usize,
    marked_bond: BondEdge,
    collapse_plain_explicit_hydrogens: bool,
) -> bool {
    let Some(center_state) = rdkit_atom_pair_conjugation_candidate_state(
        smiles,
        center_id,
        collapse_plain_explicit_hydrogens,
    ) else {
        return false;
    };
    if !(2..=3).contains(&center_state.total_degree) {
        return false;
    }

    let Some(marked_neighbor_state) = rdkit_atom_pair_conjugation_candidate_state(
        smiles,
        marked_neighbor_id,
        collapse_plain_explicit_hydrogens,
    ) else {
        return false;
    };

    let marked_is_multiple = rdkit_bond_order(marked_bond.bond_type()) >= 2;
    let marked_can_be_single_side = marked_neighbor_state.total_degree <= 3;

    smiles.edges_for_node(center_id).any(|other_edge| {
        let other_id = rdkit_atom_pair_neighbor_id(other_edge, center_id);
        if other_id == marked_neighbor_id {
            return false;
        }
        if collapse_plain_explicit_hydrogens
            && rdkit_rdk_atom_is_collapsible_hydrogen(smiles, other_id)
        {
            return false;
        }

        let Some(other_state) = rdkit_atom_pair_conjugation_candidate_state(
            smiles,
            other_id,
            collapse_plain_explicit_hydrogens,
        ) else {
            return false;
        };

        let other_is_multiple = rdkit_bond_order(other_edge.bond_type()) >= 2;
        (marked_is_multiple && other_state.total_degree <= 3)
            || (marked_can_be_single_side && other_is_multiple)
    })
}

#[inline]
fn rdkit_atom_pair_bond_is_conjugated(
    smiles: &Smiles,
    edge: BondEdge,
    collapse_plain_explicit_hydrogens: bool,
) -> bool {
    if edge.is_aromatic() {
        return true;
    }

    rdkit_atom_pair_conjugated_by_center(
        smiles,
        edge.source(),
        edge.target(),
        edge,
        collapse_plain_explicit_hydrogens,
    ) || rdkit_atom_pair_conjugated_by_center(
        smiles,
        edge.target(),
        edge.source(),
        edge,
        collapse_plain_explicit_hydrogens,
    )
}

#[inline]
fn rdkit_atom_pair_atom_has_conjugated_bond(
    smiles: &Smiles,
    atom_id: usize,
    collapse_plain_explicit_hydrogens: bool,
) -> bool {
    smiles.edges_for_node(atom_id).any(|edge| {
        let other_id = rdkit_atom_pair_neighbor_id(edge, atom_id);
        if collapse_plain_explicit_hydrogens
            && rdkit_rdk_atom_is_collapsible_hydrogen(smiles, other_id)
        {
            return false;
        }

        rdkit_atom_pair_bond_is_conjugated(smiles, edge, collapse_plain_explicit_hydrogens)
    })
}

#[inline]
fn rdkit_atom_pair_is_sp3(
    smiles: &Smiles,
    atom_id: usize,
    element: elements_rs::Element,
    state: RdkitAtomPairValenceState,
    collapse_plain_explicit_hydrogens: bool,
) -> bool {
    let orbital_count = if u8::from(element) < 89 {
        rdkit_atom_pair_num_bonds_plus_lone_pairs(element, state)
    } else {
        i32::try_from(state.total_degree).unwrap_or(i32::MAX)
    };

    orbital_count == 4
        && (state.total_degree > 3
            || !rdkit_atom_pair_atom_has_conjugated_bond(
                smiles,
                atom_id,
                collapse_plain_explicit_hydrogens,
            ))
}

#[inline]
fn rdkit_atom_pair_pi_electrons(
    smiles: &Smiles,
    atom_id: usize,
    atom: &Atom,
    collapse_plain_explicit_hydrogens: bool,
) -> u32 {
    if atom.aromatic() {
        return 1;
    }

    let Some(element) = atom.element() else {
        return 0;
    };

    if matches!(
        element,
        elements_rs::Element::Cl | elements_rs::Element::Br | elements_rs::Element::I
    ) && (hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, atom_id).is_some()
        || hydrogenated_hypervalent_halogen_center_double_oxygen_count(smiles, atom_id).is_some())
    {
        return 0;
    }

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
            ) || smiles.edge_count_for_node(other_id) < 2
            {
                continue;
            }
            if rdkit_bond_order(edge.bond_type()) >= 2
                && (hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, other_id)
                    .is_some()
                    || hydrogenated_hypervalent_halogen_center_double_oxygen_count(
                        smiles, other_id,
                    )
                    .is_some())
            {
                return 0;
            }
        }
    }

    let state =
        rdkit_atom_pair_valence_state(smiles, atom_id, atom, collapse_plain_explicit_hydrogens);
    if state.explicit_valence <= state.physical_bonds {
        return 0;
    }

    if rdkit_atom_pair_is_sp3(
        smiles,
        atom_id,
        element,
        state,
        collapse_plain_explicit_hydrogens,
    ) {
        return 0;
    }

    state.explicit_valence - state.physical_bonds
}

#[inline]
fn rdkit_atom_invariant_fields(
    smiles: &Smiles,
    atom_id: usize,
    ring_atom_membership: &RingAtomMembership,
) -> AtomInvariantFields {
    let Some(atom) = smiles.node_by_id(atom_id) else {
        return AtomInvariantFields::default();
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
    let formal_charge =
        i32::from(normalized_rdkit_charge(smiles, atom_id, atom).unwrap_or(atom.charge_value()));
    let delta_mass = morgan_delta_mass(atom);

    AtomInvariantFields {
        atomic_number,
        total_degree,
        total_hydrogens,
        formal_charge,
        delta_mass,
        in_ring: ring_atom_membership.contains_atom(atom_id),
    }
}

#[inline]
fn morgan_delta_mass(atom: &Atom) -> i32 {
    let Some(element) = atom.element() else {
        return 0;
    };

    let standard_atomic_mass = rdkit_atomic_weight(element);
    let atom_mass = morgan_isotope_mass(atom).unwrap_or(standard_atomic_mass);

    (atom_mass - standard_atomic_mass) as i32
}

#[inline]
fn morgan_isotope_mass(atom: &Atom) -> Option<f64> {
    atom.isotope_mass_number()?;
    atom.isotope()
        .ok()
        .map(|isotope| isotope.relative_atomic_mass())
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

#[derive(Debug)]
struct RdkitNormalization {
    smiles: Smiles,
}

fn normalize_smiles_for_rdkit(
    smiles: &Smiles,
) -> Result<RdkitNormalization, SmilesPreparationError> {
    let has_aromatic_labels = smiles.nodes().iter().any(Atom::aromatic)
        || (0..smiles.atom_count())
            .any(|atom_id| smiles.edges_for_node(atom_id).any(BondEdge::is_aromatic));

    let base = if has_aromatic_labels {
        if let Ok(kekule) = smiles.kekulize_standalone() {
            kekule
        } else {
            return Ok(RdkitNormalization {
                smiles: smiles.clone(),
            });
        }
    } else {
        smiles.clone()
    };

    let perception = base
        .perceive_aromaticity_for(AromaticityPolicy::RdkitDefault)
        .map_err(SmilesPreparationError::Aromaticity)?;
    Ok(RdkitNormalization {
        smiles: perception.into_aromaticized(),
    })
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

    /// Prepares a graph while preserving ordinary explicit hydrogens as
    /// fingerprint centers.
    ///
    /// This matches RDKit `AddHs`-style workflows. The default
    /// [`prepare`](Self::prepare) path mirrors sanitized `MolFromSmiles`
    /// behavior and collapses ordinary explicit hydrogens.
    #[inline]
    #[must_use]
    pub fn prepare_preserving_explicit_hydrogens<'a>(
        &'a mut self,
        smiles: &Smiles,
    ) -> SmilesRdkitGraph<'a> {
        self.try_prepare_preserving_explicit_hydrogens(smiles)
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
        self.try_prepare_with_options(smiles, true)
    }

    /// Prepares a graph while preserving ordinary explicit hydrogens as
    /// fingerprint centers.
    pub fn try_prepare_preserving_explicit_hydrogens<'a>(
        &'a mut self,
        smiles: &Smiles,
    ) -> Result<SmilesRdkitGraph<'a>, SmilesPreparationError> {
        self.try_prepare_with_options(smiles, false)
    }

    fn try_prepare_with_options<'a>(
        &'a mut self,
        smiles: &Smiles,
        collapse_plain_explicit_hydrogens: bool,
    ) -> Result<SmilesRdkitGraph<'a>, SmilesPreparationError> {
        let normalized = normalize_smiles_for_rdkit(smiles)?;
        self.normalized_smiles = Some(normalized.smiles);
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
            collapse_plain_explicit_hydrogens,
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
    collapse_plain_explicit_hydrogens: bool,
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
    fn ecfp_atom_is_ignored(&self, atom_id: usize) -> bool {
        self.collapse_plain_explicit_hydrogens
            && rdkit_rdk_atom_is_collapsible_hydrogen(self.smiles, atom_id)
    }

    #[inline]
    fn ecfp_atom_invariant_fields(&self, atom_id: usize) -> AtomInvariantFields {
        rdkit_atom_invariant_fields(self.smiles, atom_id, self.ring_atom_membership)
    }

    #[inline]
    fn ecfp_bond_invariant(&self, bond: &Self::Bond, use_bond_types: bool) -> u32 {
        rdkit_bond_code(self.smiles, *bond, use_bond_types)
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
        rdkit_atom_pair_atom_code_with_options(
            self.smiles,
            atom_id,
            self.collapse_plain_explicit_hydrogens,
        )
    }

    fn atom_pair_adjacency_and_codes(&self) -> (Vec<Vec<usize>>, Vec<u32>) {
        if !self.collapse_plain_explicit_hydrogens {
            let atom_count = self.smiles.atom_count();
            let mut adjacency = Vec::with_capacity(atom_count);
            let mut atom_codes = Vec::with_capacity(atom_count);

            for atom_id in 0..atom_count {
                let neighbors = self
                    .smiles
                    .edges_for_node(atom_id)
                    .map(|edge| {
                        if edge.source() == atom_id {
                            edge.target()
                        } else {
                            edge.source()
                        }
                    })
                    .collect();
                adjacency.push(neighbors);
                atom_codes.push(rdkit_atom_pair_atom_code(self.smiles, atom_id));
            }

            return (adjacency, atom_codes);
        }

        let atom_count = self.smiles.atom_count();
        let mut retained_atom_ids = vec![None; atom_count];
        let mut atom_codes = Vec::with_capacity(atom_count);

        for (atom_id, retained_atom_id) in retained_atom_ids.iter_mut().enumerate() {
            if rdkit_rdk_atom_is_collapsible_hydrogen(self.smiles, atom_id) {
                continue;
            }

            *retained_atom_id = Some(atom_codes.len());
            atom_codes.push(rdkit_atom_pair_atom_code_with_options(
                self.smiles,
                atom_id,
                true,
            ));
        }

        let mut adjacency = vec![Vec::new(); atom_codes.len()];
        for atom_id in 0..atom_count {
            for edge in self.smiles.edges_for_node(atom_id) {
                let other_id = if edge.source() == atom_id {
                    edge.target()
                } else {
                    edge.source()
                };
                if atom_id >= other_id {
                    continue;
                }

                let Some(source) = retained_atom_ids[atom_id] else {
                    continue;
                };
                let Some(target) = retained_atom_ids[other_id] else {
                    continue;
                };

                adjacency[source].push(target);
                adjacency[target].push(source);
            }
        }

        (adjacency, atom_codes)
    }
}

impl TopologicalTorsionGraph for Smiles {
    #[inline]
    fn topological_torsion_atom_is_hydrogen(&self, atom_id: usize) -> bool {
        rdkit_rdk_atom_is_hydrogen(self, atom_id)
    }

    #[inline]
    fn topological_torsion_atom_code(&self, atom_id: usize, branch_subtract: u8) -> u32 {
        rdkit_topological_torsion_atom_code(self, atom_id, branch_subtract)
    }
}

impl TopologicalTorsionGraph for SmilesRdkitGraph<'_> {
    #[inline]
    fn topological_torsion_atom_is_hydrogen(&self, atom_id: usize) -> bool {
        rdkit_rdk_atom_is_hydrogen(self.smiles, atom_id)
    }

    #[inline]
    fn topological_torsion_atom_code(&self, atom_id: usize, branch_subtract: u8) -> u32 {
        rdkit_topological_torsion_atom_code_with_options(
            self.smiles,
            atom_id,
            branch_subtract,
            self.collapse_plain_explicit_hydrogens,
        )
    }
}

#[inline]
fn rdkit_rdk_atom_invariant(smiles: &Smiles, atom_id: usize) -> u32 {
    let Some(atom) = smiles.node_by_id(atom_id) else {
        return 0;
    };
    let atomic_number = atom
        .element()
        .map_or(0_u32, |element| u32::from(u8::from(element)));
    ((atomic_number % 128) << 1) | u32::from(atom.aromatic())
}

#[inline]
fn rdkit_rdk_atom_is_hydrogen(smiles: &Smiles, atom_id: usize) -> bool {
    smiles.node_by_id(atom_id).and_then(Atom::element) == Some(elements_rs::Element::H)
}

#[inline]
fn rdkit_rdk_atom_is_collapsible_hydrogen(smiles: &Smiles, atom_id: usize) -> bool {
    let Some(atom) = smiles.node_by_id(atom_id) else {
        return false;
    };
    if atom.element() != Some(elements_rs::Element::H)
        || atom.isotope_mass_number().is_some()
        || atom.charge_value() != 0
    {
        return false;
    }

    let mut neighbor_ids = smiles.edges_for_node(atom_id).map(|edge| {
        if edge.source() == atom_id {
            edge.target()
        } else {
            edge.source()
        }
    });
    let Some(neighbor_id) = neighbor_ids.next() else {
        return false;
    };
    if neighbor_ids.next().is_some() {
        return false;
    }

    smiles
        .node_by_id(neighbor_id)
        .and_then(Atom::element)
        .is_some_and(|element| element != elements_rs::Element::H)
}

#[inline]
fn is_hypervalent_halogen_oxygen_bond(smiles: &Smiles, bond: BondEdge) -> bool {
    let begin_atom_id = bond.source();
    let end_atom_id = bond.target();
    let begin_atom = smiles.node_by_id(begin_atom_id);
    let end_atom = smiles.node_by_id(end_atom_id);
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

    is_halogen_oxygen_bond
        && (hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, begin_atom_id).is_some()
            || hydrogenated_hypervalent_halogen_center_double_oxygen_count(smiles, begin_atom_id)
                .is_some()
            || hypervalent_halogen_oxoacid_center_double_oxygen_count(smiles, end_atom_id)
                .is_some()
            || hydrogenated_hypervalent_halogen_center_double_oxygen_count(smiles, end_atom_id)
                .is_some())
}

#[inline]
fn is_charge_separated_phosphorylimide_oxygen_bond(smiles: &Smiles, bond: BondEdge) -> bool {
    let begin_atom_id = bond.source();
    let end_atom_id = bond.target();
    let begin_atom = smiles.node_by_id(begin_atom_id);
    let end_atom = smiles.node_by_id(end_atom_id);
    let is_phosphorus_oxygen_bond = matches!(
        (
            begin_atom.and_then(Atom::element),
            end_atom.and_then(Atom::element)
        ),
        (Some(elements_rs::Element::O), Some(elements_rs::Element::P))
            | (Some(elements_rs::Element::P), Some(elements_rs::Element::O))
    );

    is_phosphorus_oxygen_bond
        && rdkit_bond_order(bond.bond_type()) >= 2
        && (charge_separated_phosphoryl_center_double_oxygen_count(smiles, begin_atom_id).is_some()
            || charge_separated_phosphoryl_center_double_oxygen_count(smiles, end_atom_id)
                .is_some())
}

#[inline]
fn rdkit_bond_code(smiles: &Smiles, bond: BondEdge, use_bond_order: bool) -> u32 {
    if !use_bond_order
        || is_hypervalent_halogen_oxygen_bond(smiles, bond)
        || is_charge_separated_phosphorylimide_oxygen_bond(smiles, bond)
    {
        return 1;
    }

    rdkit_morgan_bond_type_code(bond)
}

#[inline]
fn rdkit_morgan_bond_type_code(bond: BondEdge) -> u32 {
    let bond_type = bond.bond_type();
    if bond_type == Bond::Triple {
        3
    } else if bond.is_aromatic() {
        12
    } else {
        rdkit_bond_order(bond_type)
    }
}

impl RdkFingerprintGraph for Smiles {
    #[inline]
    fn rdk_atom_invariant(&self, atom_id: usize) -> u32 {
        rdkit_rdk_atom_invariant(self, atom_id)
    }

    #[inline]
    fn rdk_atom_is_hydrogen(&self, atom_id: usize) -> bool {
        rdkit_rdk_atom_is_hydrogen(self, atom_id)
    }

    #[inline]
    fn rdk_atom_is_collapsible_hydrogen(&self, atom_id: usize) -> bool {
        rdkit_rdk_atom_is_collapsible_hydrogen(self, atom_id)
    }

    #[inline]
    fn rdk_bond_code(&self, bond: &Self::Bond, use_bond_order: bool) -> u32 {
        rdkit_bond_code(self, *bond, use_bond_order)
    }
}

impl RdkFingerprintGraph for SmilesRdkitGraph<'_> {
    #[inline]
    fn rdk_atom_invariant(&self, atom_id: usize) -> u32 {
        rdkit_rdk_atom_invariant(self.smiles, atom_id)
    }

    #[inline]
    fn rdk_atom_is_hydrogen(&self, atom_id: usize) -> bool {
        rdkit_rdk_atom_is_hydrogen(self.smiles, atom_id)
    }

    #[inline]
    fn rdk_atom_is_collapsible_hydrogen(&self, atom_id: usize) -> bool {
        rdkit_rdk_atom_is_collapsible_hydrogen(self.smiles, atom_id)
    }

    #[inline]
    fn rdk_bond_code(&self, bond: &Self::Bond, use_bond_order: bool) -> u32 {
        rdkit_bond_code(self.smiles, *bond, use_bond_order)
    }
}

#[inline]
fn rdkit_bond_order(bond: Bond) -> u32 {
    match bond {
        Bond::Single | Bond::Up | Bond::Down => 1,
        Bond::Double => 2,
        Bond::Triple => 3,
        Bond::Quadruple => 4,
    }
}

#[cfg(test)]
mod tests {
    use alloc::{vec, vec::Vec};

    use elements_rs::Element;
    use geometric_traits::traits::{Edges, Graph, MonoplexGraph, TypedNode};
    use smiles_parser::{
        atom::Atom,
        bond::{Bond, bond_edge::BondEdge},
        smiles::{AromaticityAssignmentApplicationError, Smiles},
    };

    use super::{
        SmilesPreparationError, SmilesRdkitScratch,
        charge_separated_phosphoryl_center_double_oxygen_count,
        hypervalent_halogen_oxoacid_center_double_oxygen_count, morgan_delta_mass,
        normalize_smiles_for_rdkit, normalized_rdkit_charge, rdkit_atom_invariant_fields,
        rdkit_atom_pair_atom_code, rdkit_atom_pair_pi_electrons, rdkit_atom_pair_type_index,
        rdkit_atomic_weight, rdkit_bond_order, rdkit_morgan_bond_type_code,
        rdkit_topological_torsion_atom_code,
    };
    use crate::{AtomPairGraph, EcfpGraph, MolecularAtom, MolecularGraph, TopologicalTorsionGraph};

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
            ("[3H]C", vec![4277593706, 2246733040]),
            (
                "C(=[P+]=O)C(=O)N",
                vec![
                    2246703798, 1010215739, 864942730, 2246699815, 864942730, 847957139,
                ],
            ),
            (
                "C1=CC=P(=O)C=C1",
                vec![
                    3218693969, 3218693969, 3218693969, 1660341273, 864942795, 3218693969,
                    3218693969,
                ],
            ),
            (
                "c1op[o+]c1",
                vec![3218693969, 3189457552, 2261158082, 3189554341, 3218693969],
            ),
            (
                "C1=CC=C(C=C1)[IH](=O)C2=CC=CC=C2",
                vec![
                    3218693969, 3218693969, 3218693969, 3217380708, 3218693969, 3218693969,
                    322229836, 864942730, 3217380708, 3218693969, 3218693969, 3218693969,
                    3218693969, 3218693969,
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
    fn connectivity_invariants_use_elements_rs_masses_for_known_rdkit_isotope_errors() {
        for (smiles, expected) in [("[162Tm]", vec![3275783703]), ("[239Th]", vec![496497324])] {
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

    #[test]
    fn smiles_parser_adapter_helpers_cover_atom_and_bond_views() {
        let smiles: Smiles = "C=O".parse().expect("fixture SMILES should parse");
        let atom = smiles.node_by_id(0).expect("first atom should exist");
        let bond = smiles
            .edge_for_node_pair((0, 1))
            .expect("fixture bond should exist");

        assert_eq!(atom.atom_type(), atom.node_type());
        assert_eq!(bond.source(), 0);
        assert_eq!(bond.target(), 1);
        assert_eq!(bond.bond_type(), smiles_parser::bond::Bond::Double);
        assert_eq!(MolecularGraph::atom(&smiles, 0), Some(atom));
        assert_eq!(MolecularGraph::bonds(&smiles, 0).count(), 1);
    }

    #[test]
    fn atomic_weight_and_delta_mass_helpers_match_expected_values() {
        for (element, expected) in [
            (Element::Tc, 98.0),
            (Element::Pm, 145.0),
            (Element::Po, 209.0),
            (Element::At, 210.0),
            (Element::Rn, 222.0),
            (Element::Fr, 223.0),
            (Element::Ra, 226.0),
            (Element::Ac, 227.0),
            (Element::Np, 237.0),
            (Element::Pu, 244.0),
            (Element::Am, 243.0),
            (Element::Cm, 247.0),
            (Element::Bk, 247.0),
            (Element::Cf, 251.0),
            (Element::Es, 252.0),
            (Element::Fm, 257.0),
            (Element::Md, 258.0),
            (Element::No, 259.0),
            (Element::Lr, 262.0),
        ] {
            assert_eq!(rdkit_atomic_weight(element), expected);
        }
        assert_eq!(rdkit_atomic_weight(Element::Cl), 35.45);

        let unlabeled: Smiles = "CC".parse().expect("fixture SMILES should parse");
        let isotope: Smiles = "[13CH4]".parse().expect("fixture SMILES should parse");
        let mercury: Smiles = "[Hg]".parse().expect("fixture SMILES should parse");
        let thulium_162: Smiles = "[162Tm]".parse().expect("fixture SMILES should parse");
        let thorium_239: Smiles = "[239Th]".parse().expect("fixture SMILES should parse");

        assert_eq!(morgan_delta_mass(unlabeled.node_by_id(0).unwrap()), 0);
        assert_eq!(morgan_delta_mass(isotope.node_by_id(0).unwrap()), 0);
        assert_eq!(morgan_delta_mass(mercury.node_by_id(0).unwrap()), 0);
        assert_eq!(morgan_delta_mass(thulium_162.node_by_id(0).unwrap()), -6);
        assert_eq!(morgan_delta_mass(thorium_239.node_by_id(0).unwrap()), 7);
    }

    #[test]
    fn normalized_charge_helpers_detect_rdkit_special_cases() {
        let smiles: Smiles = "O=Cl(=O)O".parse().expect("fixture SMILES should parse");
        let non_halogen: Smiles = "CO".parse().expect("fixture SMILES should parse");
        let phosphorus: Smiles = "CN=P(=O)O".parse().expect("fixture SMILES should parse");
        let phosphorus_carbon: Smiles = "CC=P(=O)C".parse().expect("fixture SMILES should parse");

        assert_eq!(
            hypervalent_halogen_oxoacid_center_double_oxygen_count(&smiles, 1),
            Some(2)
        );
        assert_eq!(
            normalized_rdkit_charge(&smiles, 1, smiles.node_by_id(1).unwrap()),
            Some(2)
        );
        assert_eq!(
            normalized_rdkit_charge(&smiles, 0, smiles.node_by_id(0).unwrap()),
            Some(-1)
        );
        assert_eq!(
            hypervalent_halogen_oxoacid_center_double_oxygen_count(&smiles, 0),
            None
        );
        assert_eq!(
            hypervalent_halogen_oxoacid_center_double_oxygen_count(&non_halogen, 0),
            None
        );

        assert_eq!(
            charge_separated_phosphoryl_center_double_oxygen_count(&phosphorus, 2),
            Some(1)
        );
        assert_eq!(
            normalized_rdkit_charge(&phosphorus, 2, phosphorus.node_by_id(2).unwrap()),
            Some(1)
        );
        assert_eq!(
            (0..phosphorus.atom_count())
                .filter(
                    |&atom_id| phosphorus.node_by_id(atom_id).and_then(Atom::element)
                        == Some(Element::O)
                )
                .filter_map(|atom_id| {
                    normalized_rdkit_charge(
                        &phosphorus,
                        atom_id,
                        phosphorus.node_by_id(atom_id).unwrap(),
                    )
                })
                .collect::<Vec<_>>(),
            vec![-1]
        );
        assert_eq!(
            charge_separated_phosphoryl_center_double_oxygen_count(&phosphorus_carbon, 2),
            Some(1)
        );
        assert_eq!(
            normalized_rdkit_charge(
                &phosphorus_carbon,
                2,
                phosphorus_carbon.node_by_id(2).unwrap()
            ),
            Some(1)
        );
    }

    #[test]
    fn rdkit_atom_pair_type_index_covers_known_and_wildcard_atoms() {
        let carbon: Smiles = "C".parse().expect("fixture SMILES should parse");
        let iron: Smiles = "[Fe]".parse().expect("fixture SMILES should parse");

        assert_eq!(rdkit_atom_pair_type_index(carbon.node_by_id(0).unwrap()), 1);
        assert_eq!(
            rdkit_atom_pair_type_index(iron.node_by_id(0).unwrap()),
            u32::try_from(15).unwrap()
        );
    }

    #[test]
    fn rdkit_atom_code_helpers_cover_invalid_atoms_and_special_pi_cases() {
        let aromatic: Smiles = "c1ccccc1".parse().expect("fixture SMILES should parse");
        let halogen: Smiles = "O=Cl(=O)O".parse().expect("fixture SMILES should parse");
        let bromine: Smiles = "O=Br(=O)O".parse().expect("fixture SMILES should parse");
        let iodine: Smiles = "O=I(=O)O".parse().expect("fixture SMILES should parse");
        let phosphorus: Smiles = "OP(=O)O".parse().expect("fixture SMILES should parse");
        let sulfur: Smiles = "O=S=O".parse().expect("fixture SMILES should parse");
        let tungsten: Smiles = "O=[W](=O)(O)O"
            .parse()
            .expect("fixture SMILES should parse");

        assert_eq!(rdkit_atom_pair_atom_code(&aromatic, 99), 0);
        assert_eq!(rdkit_topological_torsion_atom_code(&aromatic, 99, 1), 0);
        assert_eq!(
            rdkit_atom_pair_pi_electrons(&aromatic, 0, aromatic.node_by_id(0).unwrap(), false),
            1
        );
        assert_eq!(
            rdkit_atom_pair_pi_electrons(&halogen, 0, halogen.node_by_id(0).unwrap(), false),
            0
        );
        assert_eq!(
            rdkit_atom_pair_pi_electrons(&halogen, 1, halogen.node_by_id(1).unwrap(), false),
            0
        );
        assert_eq!(
            rdkit_atom_pair_pi_electrons(&bromine, 1, bromine.node_by_id(1).unwrap(), false),
            0
        );
        assert_eq!(
            rdkit_atom_pair_pi_electrons(&iodine, 1, iodine.node_by_id(1).unwrap(), false),
            0
        );
        assert_eq!(
            rdkit_atom_pair_pi_electrons(&phosphorus, 1, phosphorus.node_by_id(1).unwrap(), false),
            0
        );
        assert_eq!(
            rdkit_atom_pair_pi_electrons(&sulfur, 1, sulfur.node_by_id(1).unwrap(), false),
            2
        );
        assert_eq!(
            rdkit_atom_pair_pi_electrons(&tungsten, 1, tungsten.node_by_id(1).unwrap(), false),
            0
        );
    }

    #[test]
    fn rdkit_atom_pair_metal_atom_codes_match_reference_fixtures() {
        let expected = vec![
            ("C=CC=[Fe]", 3, 481, 480),
            ("CCOC=[Fe]", 4, 481, 480),
            ("C=C=C=[Rh]", 3, 481, 480),
            ("O=[Sm]O[Sm]=O", 1, 490, 489),
            ("O=[Eu]O[Eu]=O", 1, 482, 481),
        ];
        let observed = expected
            .iter()
            .map(|(smiles, atom_id, _, _)| {
                let original_smiles = *smiles;
                let smiles: Smiles = original_smiles
                    .parse()
                    .expect("fixture SMILES should parse");
                let mut scratch = SmilesRdkitScratch::default();
                let graph = scratch.prepare(&smiles);

                (
                    original_smiles,
                    *atom_id,
                    AtomPairGraph::atom_pair_atom_code(&graph, *atom_id),
                    TopologicalTorsionGraph::topological_torsion_atom_code(&graph, *atom_id, 1),
                )
            })
            .collect::<Vec<_>>();

        assert_eq!(observed, expected);
    }

    #[test]
    fn rdkit_atom_pair_metal_pi_electron_counts_match_reference_fixtures() {
        let expected = vec![
            ("C=CC=[Fe]", 3, 0),
            ("CCOC=[Fe]", 4, 0),
            ("C=C=C=[Rh]", 3, 0),
            ("O=[Sm]O[Sm]=O", 1, 1),
            ("O=[Eu]O[Eu]=O", 1, 0),
        ];
        let observed = expected
            .iter()
            .map(|(smiles, atom_id, _)| {
                let original_smiles = *smiles;
                let smiles: Smiles = original_smiles
                    .parse()
                    .expect("fixture SMILES should parse");
                let mut scratch = SmilesRdkitScratch::default();
                let graph = scratch.prepare(&smiles);
                let atom =
                    MolecularGraph::atom(&graph, *atom_id).expect("fixture atom should exist");

                (
                    original_smiles,
                    *atom_id,
                    rdkit_atom_pair_pi_electrons(graph.smiles, *atom_id, atom, true),
                )
            })
            .collect::<Vec<_>>();

        assert_eq!(observed, expected);
    }

    #[test]
    fn rdkit_atom_invariant_fields_handles_invalid_atom_ids() {
        // Out-of-range atom ids must not panic; the field builder returns
        // the default tuple (all zeros, in_ring = false) so downstream code
        // sees a deterministic, well-defined value.
        let smiles: Smiles = "CC".parse().expect("fixture SMILES should parse");
        let ring = smiles.ring_atom_membership();
        assert_eq!(
            rdkit_atom_invariant_fields(&smiles, 99, &ring),
            crate::AtomInvariantFields::default(),
        );
    }

    #[test]
    fn smiles_preparation_error_display_includes_inner_error() {
        let error = SmilesPreparationError::Aromaticity(
            AromaticityAssignmentApplicationError::MissingBondEdge {
                node_a: 1,
                node_b: 2,
            },
        );

        let rendered = alloc::format!("{error}");
        assert!(rendered.contains("failed to apply RDKit-default aromaticity"));
        assert!(rendered.contains("bond edge does not exist in graph"));
    }

    #[test]
    fn normalize_and_prepare_smiles_for_rdkit_cover_graph_wrapper_methods() {
        let aromatic: Smiles = "c1ccccc1".parse().expect("fixture SMILES should parse");
        let normalized =
            normalize_smiles_for_rdkit(&aromatic).expect("normalization should succeed");
        assert!(normalized.smiles.nodes().iter().all(|atom| atom.aromatic()));

        let mut scratch = SmilesRdkitScratch::default();
        let smiles: Smiles = "C=O".parse().expect("fixture SMILES should parse");
        let graph = scratch
            .try_prepare(&smiles)
            .expect("preparation should succeed");

        assert!(graph.has_nodes());
        assert!(graph.has_edges());
        assert_eq!(Edges::number_of_edges(graph.edges()), 2);
        assert_eq!(
            MolecularGraph::atom(&graph, 0).and_then(|atom| atom.element()),
            Some(Element::C)
        );
        assert_eq!(MolecularGraph::bonds(&graph, 0).count(), 1);
        assert_ne!(
            TopologicalTorsionGraph::topological_torsion_atom_code(&smiles, 0, 1),
            0
        );
        assert_ne!(AtomPairGraph::atom_pair_atom_code(&graph, 0), 0);
        assert_ne!(
            TopologicalTorsionGraph::topological_torsion_atom_code(&graph, 0, 1),
            0
        );

        let double_bond = MolecularGraph::bonds(&graph, 0)
            .next()
            .expect("fixture bond should exist");
        assert_eq!(graph.ecfp_bond_invariant(&double_bond, false), 1);
        assert_eq!(graph.ecfp_bond_invariant(&double_bond, true), 2);

        let halogen_smiles: Smiles = "O=Cl(=O)O".parse().expect("fixture SMILES should parse");
        let mut halogen_scratch = SmilesRdkitScratch::default();
        let halogen_graph = halogen_scratch.prepare(&halogen_smiles);
        let halogen_bond = BondEdge::new(0, 1, Bond::Double, None);
        assert_eq!(halogen_graph.ecfp_bond_invariant(&halogen_bond, true), 1);

        let quadruple_bond = BondEdge::new(0, 1, Bond::Quadruple, None);
        assert_eq!(graph.ecfp_bond_invariant(&quadruple_bond, true), 4);
    }

    #[test]
    fn rdkit_bond_order_covers_all_supported_bond_variants() {
        use smiles_parser::bond::Bond;

        assert_eq!(rdkit_bond_order(Bond::Single), 1);
        assert_eq!(rdkit_bond_order(Bond::Up), 1);
        assert_eq!(rdkit_bond_order(Bond::Down), 1);
        assert_eq!(rdkit_bond_order(Bond::Double), 2);
        assert_eq!(rdkit_bond_order(Bond::Triple), 3);
        assert_eq!(rdkit_bond_order(Bond::Quadruple), 4);
    }

    #[test]
    fn rdkit_morgan_bond_type_code_prefers_triple_over_aromaticity() {
        assert_eq!(
            rdkit_morgan_bond_type_code(BondEdge::with_aromaticity(0, 1, Bond::Triple, None, true)),
            3
        );
        assert_eq!(
            rdkit_morgan_bond_type_code(BondEdge::with_aromaticity(0, 1, Bond::Single, None, true)),
            12
        );
        assert_eq!(
            rdkit_morgan_bond_type_code(BondEdge::new(0, 1, Bond::Double, None)),
            2
        );
    }
}
