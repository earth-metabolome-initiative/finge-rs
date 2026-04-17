use alloc::vec::Vec;

use geometric_traits::traits::{Graph, MonopartiteGraph, MonoplexGraph, Vocabulary};

/// Trait implemented by molecular atom node types.
pub trait MolecularAtom {
    /// Chemistry-facing atom identity used by simple fingerprints.
    type AtomType: Clone + Eq;

    /// Returns the chemistry-facing atom identity for this node.
    fn atom_type(&self) -> Self::AtomType;
}

/// Trait implemented by molecular bond edge values.
pub trait MolecularBond {
    /// Dense node identifier type used by the parent graph.
    type NodeId: Copy + Eq;
    /// Chemistry-facing bond identity.
    type BondType: Clone + Eq;

    /// Returns the source atom identifier.
    fn source(&self) -> Self::NodeId;

    /// Returns the destination atom identifier.
    fn target(&self) -> Self::NodeId;

    /// Returns the chemistry-facing bond identity.
    fn bond_type(&self) -> Self::BondType;
}

/// Trait implemented by molecular graphs backed by `geometric-traits`.
pub trait MolecularGraph: Graph + MonopartiteGraph + MonoplexGraph
where
    Self::NodeSymbol: MolecularAtom,
{
    /// Bond view returned when traversing incident edges.
    type Bond: MolecularBond<NodeId = Self::NodeId>;

    /// Returns the atom for a node identifier, if present.
    fn atom(&self, node_id: Self::NodeId) -> Option<&Self::NodeSymbol>;

    /// Returns the incident bonds for a node identifier.
    fn bonds(&self, node_id: Self::NodeId) -> impl Iterator<Item = Self::Bond> + '_;

    /// Returns the number of atoms in the graph.
    #[inline]
    fn atom_count(&self) -> usize {
        self.nodes_vocabulary().len()
    }

    /// Returns whether the graph has no atoms.
    #[inline]
    fn is_empty_molecule(&self) -> bool {
        self.atom_count() == 0
    }
}

/// Graphs that can provide RDKit-style ECFP connectivity invariants.
pub trait EcfpGraph: MolecularGraph<NodeId = usize>
where
    Self::NodeSymbol: MolecularAtom,
{
    /// Returns the RDKit-style atom invariant used to seed ECFP.
    fn ecfp_atom_invariant(&self, atom_id: usize, include_ring_membership: bool) -> u32;

    /// Returns the RDKit-style bond invariant used during neighborhood expansion.
    fn ecfp_bond_invariant(&self, bond: &Self::Bond, use_bond_types: bool) -> u32;
}

/// Graphs that can provide RDKit-style AtomPair atom codes.
pub trait AtomPairGraph: MolecularGraph<NodeId = usize>
where
    Self::NodeSymbol: MolecularAtom,
{
    /// Returns the RDKit-style AtomPair atom code for one atom.
    fn atom_pair_atom_code(&self, atom_id: usize) -> u32;

    /// Returns the per-atom adjacency lists and RDKit-style AtomPair atom
    /// codes.
    ///
    /// Backends can override this to build both in one pass when that is
    /// cheaper than computing them separately.
    fn atom_pair_adjacency_and_codes(&self) -> (Vec<Vec<usize>>, Vec<u32>)
    where
        Self::Bond: MolecularBond<NodeId = usize>,
    {
        let atom_count = self.atom_count();
        let mut adjacency = Vec::with_capacity(atom_count);
        let mut atom_codes = Vec::with_capacity(atom_count);

        for atom_id in 0..atom_count {
            let neighbors = self
                .bonds(atom_id)
                .filter_map(|bond| {
                    if bond.source() == atom_id {
                        Some(bond.target())
                    } else if bond.target() == atom_id {
                        Some(bond.source())
                    } else {
                        None
                    }
                })
                .collect();
            adjacency.push(neighbors);
            atom_codes.push(self.atom_pair_atom_code(atom_id));
        }

        (adjacency, atom_codes)
    }
}
