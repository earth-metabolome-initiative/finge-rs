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

/// Graphs that can provide RDKit-style Topological Torsion atom codes.
pub trait TopologicalTorsionGraph: MolecularGraph<NodeId = usize>
where
    Self::NodeSymbol: MolecularAtom,
{
    /// Returns the RDKit-style Topological Torsion atom code for one atom.
    ///
    /// `branch_subtract` matches RDKit's `AtomPairs::getAtomCode()` behavior
    /// for torsions: endpoints subtract `1` and internal atoms subtract `2`.
    fn topological_torsion_atom_code(&self, atom_id: usize, branch_subtract: u8) -> u32;

    /// Returns the per-atom adjacency lists for path enumeration.
    #[inline]
    fn topological_torsion_adjacency(&self) -> Vec<Vec<usize>>
    where
        Self::Bond: MolecularBond<NodeId = usize>,
    {
        let atom_count = self.atom_count();
        let mut adjacency = Vec::with_capacity(atom_count);

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
        }

        adjacency
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec;

    use geometric_traits::traits::{Graph, MonopartiteGraph, MonoplexGraph};
    use smiles_parser::{
        bond::{Bond, bond_edge::BondEdge},
        smiles::Smiles,
    };

    use crate::traits::{AtomPairGraph, MolecularGraph, TopologicalTorsionGraph};

    struct MalformedBondSmiles {
        inner: Smiles,
        malformed_on: usize,
        malformed_bond: BondEdge,
    }

    impl MalformedBondSmiles {
        fn new(smiles: &str, malformed_on: usize, malformed_bond: BondEdge) -> Self {
            Self {
                inner: smiles.parse().expect("fixture SMILES should parse"),
                malformed_on,
                malformed_bond,
            }
        }
    }

    impl Graph for MalformedBondSmiles {
        fn has_nodes(&self) -> bool {
            self.inner.has_nodes()
        }

        fn has_edges(&self) -> bool {
            self.inner.has_edges()
        }
    }

    impl MonoplexGraph for MalformedBondSmiles {
        type Edge = <Smiles as MonoplexGraph>::Edge;
        type Edges = <Smiles as MonoplexGraph>::Edges;

        fn edges(&self) -> &Self::Edges {
            self.inner.edges()
        }
    }

    impl MonopartiteGraph for MalformedBondSmiles {
        type NodeId = <Smiles as MonopartiteGraph>::NodeId;
        type NodeSymbol = <Smiles as MonopartiteGraph>::NodeSymbol;
        type Nodes = <Smiles as MonopartiteGraph>::Nodes;

        fn nodes_vocabulary(&self) -> &Self::Nodes {
            self.inner.nodes_vocabulary()
        }
    }

    impl MolecularGraph for MalformedBondSmiles {
        type Bond = BondEdge;

        fn atom(&self, node_id: Self::NodeId) -> Option<&Self::NodeSymbol> {
            self.inner.node_by_id(node_id)
        }

        fn bonds(&self, node_id: Self::NodeId) -> impl Iterator<Item = Self::Bond> + '_ {
            let mut bonds = self
                .inner
                .edges_for_node(node_id)
                .collect::<alloc::vec::Vec<_>>();
            if node_id == self.malformed_on {
                bonds.push(self.malformed_bond);
            }
            bonds.into_iter()
        }
    }

    impl AtomPairGraph for MalformedBondSmiles {
        fn atom_pair_atom_code(&self, atom_id: usize) -> u32 {
            atom_id as u32 + 1
        }
    }

    impl TopologicalTorsionGraph for MalformedBondSmiles {
        fn topological_torsion_atom_code(&self, atom_id: usize, branch_subtract: u8) -> u32 {
            atom_id as u32 + u32::from(branch_subtract)
        }
    }

    #[test]
    fn molecular_graph_default_helpers_cover_empty_and_non_empty_smiles() {
        let empty = Smiles::new();
        assert_eq!(MolecularGraph::atom_count(&empty), 0);
        assert!(MolecularGraph::is_empty_molecule(&empty));

        let smiles: Smiles = "CCO".parse().expect("fixture SMILES should parse");
        assert_eq!(MolecularGraph::atom_count(&smiles), 3);
        assert!(!MolecularGraph::is_empty_molecule(&smiles));
    }

    #[test]
    fn atom_pair_default_adjacency_and_codes_match_smiles_topology() {
        let empty = Smiles::new();
        let (empty_adjacency, empty_codes) = AtomPairGraph::atom_pair_adjacency_and_codes(&empty);
        assert!(empty_adjacency.is_empty());
        assert!(empty_codes.is_empty());

        let smiles: Smiles = "CCO".parse().expect("fixture SMILES should parse");
        let (mut adjacency, atom_codes) = AtomPairGraph::atom_pair_adjacency_and_codes(&smiles);
        for neighbors in &mut adjacency {
            neighbors.sort_unstable();
        }

        let expected_codes = (0..MolecularGraph::atom_count(&smiles))
            .map(|atom_id| AtomPairGraph::atom_pair_atom_code(&smiles, atom_id))
            .collect::<alloc::vec::Vec<_>>();

        assert_eq!(adjacency, vec![vec![1], vec![0, 2], vec![1]]);
        assert_eq!(atom_codes, expected_codes);
    }

    #[test]
    fn topological_torsion_default_adjacency_matches_smiles_topology() {
        let empty = Smiles::new();
        assert!(TopologicalTorsionGraph::topological_torsion_adjacency(&empty).is_empty());

        let smiles: Smiles = "C1CC1O".parse().expect("fixture SMILES should parse");
        let mut adjacency = TopologicalTorsionGraph::topological_torsion_adjacency(&smiles);
        for neighbors in &mut adjacency {
            neighbors.sort_unstable();
        }

        assert_eq!(
            adjacency,
            vec![vec![1, 2], vec![0, 2], vec![0, 1, 3], vec![2]]
        );
    }

    #[test]
    fn default_adjacency_helpers_ignore_malformed_bonds() {
        let malformed_bond: BondEdge = (1, 2, Bond::Single, None);
        let graph = MalformedBondSmiles::new("CCO", 0, malformed_bond);

        let (mut atom_pair_adjacency, atom_codes) =
            AtomPairGraph::atom_pair_adjacency_and_codes(&graph);
        for neighbors in &mut atom_pair_adjacency {
            neighbors.sort_unstable();
        }
        assert_eq!(atom_pair_adjacency, vec![vec![1], vec![0, 2], vec![1]]);
        assert_eq!(atom_codes, vec![1, 2, 3]);

        let mut torsion_adjacency = TopologicalTorsionGraph::topological_torsion_adjacency(&graph);
        for neighbors in &mut torsion_adjacency {
            neighbors.sort_unstable();
        }
        assert_eq!(torsion_adjacency, vec![vec![1], vec![0, 2], vec![1]]);
    }
}
