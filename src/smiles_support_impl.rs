use elements_rs::isotopes::RelativeAtomicMass;
use geometric_traits::traits::{Graph, MonopartiteGraph, MonoplexGraph, TypedNode};
use smiles_parser::{
    atom::{Atom, McesAtomType},
    bond::{Bond, bond_edge::BondEdge},
    smiles::{RingAtomMembership, RingAtomMembershipScratch, Smiles},
};

use crate::traits::{EcfpGraph, MolecularAtom, MolecularBond, MolecularGraph};

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
fn rdkit_atom_invariant(
    smiles: &Smiles,
    atom_id: usize,
    include_ring_membership: bool,
    ring_atom_membership: &RingAtomMembership,
) -> u32 {
    let Some(atom) = smiles.node_by_id(atom_id) else {
        return 0;
    };

    let explicit_hydrogens = u32::from(atom.hydrogen_count());
    let implicit_hydrogens = u32::from(smiles.implicit_hydrogen_count(atom_id));
    let total_hydrogens = explicit_hydrogens + implicit_hydrogens;
    let total_degree = smiles.out_degree(atom_id) as u32 + total_hydrogens;
    let atomic_number = atom
        .element()
        .map_or(0, |element| u32::from(u8::from(element)));
    let formal_charge = atom.charge_value() as u32;
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

    let standard_atomic_weight = element.standard_atomic_weight();
    let atom_mass = atom
        .isotope()
        .map(|isotope| isotope.relative_atomic_mass())
        .unwrap_or(standard_atomic_weight);

    (atom_mass - standard_atomic_weight) as i32
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

/// Reusable scratch for batch ECFP processing over `smiles-parser` graphs.
///
/// This keeps the atom-only ring-membership output and DFS scratch buffers
/// alive across molecules so callers can avoid reallocating them on every
/// fingerprint call.
#[derive(Debug, Default)]
pub struct SmilesEcfpScratch {
    ring_atom_membership: RingAtomMembership,
    ring_atom_membership_scratch: RingAtomMembershipScratch,
}

impl SmilesEcfpScratch {
    /// Prepares a `smiles-parser` graph for ECFP computation while reusing the
    /// atom-ring-membership buffers held by this scratch object.
    #[inline]
    #[must_use]
    pub fn prepare<'a>(&'a mut self, smiles: &'a Smiles) -> SmilesEcfpGraph<'a> {
        smiles.write_ring_atom_membership(
            &mut self.ring_atom_membership,
            &mut self.ring_atom_membership_scratch,
        );

        SmilesEcfpGraph {
            smiles,
            ring_atom_membership: &self.ring_atom_membership,
        }
    }
}

/// Scratch-prepared `smiles-parser` graph view for ECFP computation.
///
/// Create this through [`SmilesEcfpScratch::prepare`] and pass it to the
/// normal [`Fingerprint::compute`](crate::Fingerprint::compute) path.
#[derive(Debug, Clone, Copy)]
pub struct SmilesEcfpGraph<'a> {
    smiles: &'a Smiles,
    ring_atom_membership: &'a RingAtomMembership,
}

impl Graph for SmilesEcfpGraph<'_> {
    #[inline]
    fn has_nodes(&self) -> bool {
        self.smiles.has_nodes()
    }

    #[inline]
    fn has_edges(&self) -> bool {
        self.smiles.has_edges()
    }
}

impl MonoplexGraph for SmilesEcfpGraph<'_> {
    type Edge = <Smiles as MonoplexGraph>::Edge;
    type Edges = <Smiles as MonoplexGraph>::Edges;

    #[inline]
    fn edges(&self) -> &Self::Edges {
        self.smiles.edges()
    }
}

impl MonopartiteGraph for SmilesEcfpGraph<'_> {
    type NodeId = <Smiles as MonopartiteGraph>::NodeId;
    type NodeSymbol = <Smiles as MonopartiteGraph>::NodeSymbol;
    type Nodes = <Smiles as MonopartiteGraph>::Nodes;

    #[inline]
    fn nodes_vocabulary(&self) -> &Self::Nodes {
        self.smiles.nodes_vocabulary()
    }
}

impl MolecularGraph for SmilesEcfpGraph<'_> {
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

impl EcfpGraph for SmilesEcfpGraph<'_> {
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

        match bond.bond_type() {
            Bond::Single | Bond::Up | Bond::Down => 1,
            Bond::Double => 2,
            Bond::Triple => 3,
            Bond::Quadruple => 4,
            Bond::Aromatic => 12,
        }
    }
}

#[cfg(test)]
mod tests {
    use alloc::{vec, vec::Vec};

    use smiles_parser::smiles::Smiles;

    use super::SmilesEcfpScratch;
    use crate::{EcfpGraph, MolecularGraph};

    fn observed_connectivity_invariants(smiles: &str, include_ring_membership: bool) -> Vec<u32> {
        let smiles: Smiles = smiles.parse().expect("fixture SMILES should parse");
        let mut scratch = SmilesEcfpScratch::default();
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
