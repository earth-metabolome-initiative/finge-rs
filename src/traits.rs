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
