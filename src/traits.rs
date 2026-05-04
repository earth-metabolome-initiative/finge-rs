use alloc::{collections::BTreeMap, vec, vec::Vec};

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
    /// Returns whether this atom should be ignored as an ECFP center.
    ///
    /// Backends that mimic RDKit's sanitized SMILES path can use this to
    /// suppress ordinary explicit hydrogens that RDKit removes during parsing.
    #[inline]
    fn ecfp_atom_is_ignored(&self, _atom_id: usize) -> bool {
        false
    }

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
    /// Returns whether this atom should be excluded from torsion paths.
    ///
    /// RDKit's default Topological Torsion generator calls its path finder with
    /// `useHs=false`, so explicit hydrogen atoms are not path atoms.
    #[inline]
    fn topological_torsion_atom_is_hydrogen(&self, _atom_id: usize) -> bool {
        false
    }

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
            if self.topological_torsion_atom_is_hydrogen(atom_id) {
                adjacency.push(Vec::new());
                continue;
            }

            let neighbors = self
                .bonds(atom_id)
                .filter_map(|bond| {
                    let neighbor = if bond.source() == atom_id {
                        Some(bond.target())
                    } else if bond.target() == atom_id {
                        Some(bond.source())
                    } else {
                        None
                    }?;
                    (!self.topological_torsion_atom_is_hydrogen(neighbor)).then_some(neighbor)
                })
                .collect();
            adjacency.push(neighbors);
        }

        adjacency
    }
}

/// Bond record used by [`RdkFingerprintGraph`] path enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RdkFingerprintBond {
    /// Start atom identifier.
    pub source: usize,
    /// End atom identifier.
    pub target: usize,
    /// RDKit-style bond code used in bond hashing.
    pub bond_code: u32,
}

/// Graphs that can provide RDKit-style topological path fingerprint inputs.
pub trait RdkFingerprintGraph: MolecularGraph<NodeId = usize>
where
    Self::NodeSymbol: MolecularAtom,
{
    /// Returns the RDKit path-fingerprint atom invariant for one atom.
    fn rdk_atom_invariant(&self, atom_id: usize) -> u32;

    /// Returns whether the atom is hydrogen.
    fn rdk_atom_is_hydrogen(&self, atom_id: usize) -> bool;

    /// Returns whether the atom is an ordinary explicit hydrogen that RDKit
    /// collapses onto a neighboring heavy atom during `MolFromSmiles()`.
    #[inline]
    fn rdk_atom_is_collapsible_hydrogen(&self, _atom_id: usize) -> bool {
        false
    }

    /// Returns the RDKit path-fingerprint bond code for one bond.
    fn rdk_bond_code(&self, bond: &Self::Bond, use_bond_order: bool) -> u32;

    /// Returns atom invariants, a unique bond list, and per-atom bond
    /// adjacency for RDKFingerprint path enumeration.
    fn rdk_fingerprint_graph(
        &self,
        use_hs: bool,
        use_bond_order: bool,
    ) -> (Vec<u32>, Vec<RdkFingerprintBond>, Vec<Vec<usize>>)
    where
        Self::Bond: MolecularBond<NodeId = usize>,
    {
        let atom_count = self.atom_count();
        let atom_invariants = (0..atom_count)
            .map(|atom_id| self.rdk_atom_invariant(atom_id))
            .collect::<Vec<_>>();

        let mut unique_bonds = BTreeMap::new();
        for atom_id in 0..atom_count {
            for bond in self.bonds(atom_id) {
                let (source, target) = (bond.source(), bond.target());
                let Some(other) = (if source == atom_id {
                    Some(target)
                } else if target == atom_id {
                    Some(source)
                } else {
                    None
                }) else {
                    continue;
                };

                if self.rdk_atom_is_collapsible_hydrogen(atom_id)
                    || self.rdk_atom_is_collapsible_hydrogen(other)
                {
                    continue;
                }

                if !use_hs
                    && (self.rdk_atom_is_hydrogen(atom_id) || self.rdk_atom_is_hydrogen(other))
                {
                    continue;
                }

                let key = if atom_id <= other {
                    (atom_id, other)
                } else {
                    (other, atom_id)
                };
                unique_bonds
                    .entry(key)
                    .or_insert_with(|| self.rdk_bond_code(&bond, use_bond_order));
            }
        }

        let mut bonds = Vec::with_capacity(unique_bonds.len());
        let mut atom_bonds = vec![Vec::new(); atom_count];
        for ((source, target), bond_code) in unique_bonds {
            let bond_id = bonds.len();
            bonds.push(RdkFingerprintBond {
                source,
                target,
                bond_code,
            });
            atom_bonds[source].push(bond_id);
            atom_bonds[target].push(bond_id);
        }

        (atom_invariants, bonds, atom_bonds)
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

    use crate::traits::{
        AtomPairGraph, MolecularGraph, RdkFingerprintGraph, TopologicalTorsionGraph,
    };

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

    impl RdkFingerprintGraph for MalformedBondSmiles {
        fn rdk_atom_invariant(&self, atom_id: usize) -> u32 {
            atom_id as u32 + 10
        }

        fn rdk_atom_is_hydrogen(&self, atom_id: usize) -> bool {
            atom_id == 2
        }

        fn rdk_bond_code(&self, _bond: &Self::Bond, use_bond_order: bool) -> u32 {
            if use_bond_order { 7 } else { 1 }
        }
    }

    #[test]
    fn molecular_graph_default_helpers_cover_empty_and_non_empty_smiles() {
        let smiles: Smiles = "CCO".parse().expect("fixture SMILES should parse");
        assert_eq!(MolecularGraph::atom_count(&smiles), 3);
        assert!(!MolecularGraph::is_empty_molecule(&smiles));
    }

    #[test]
    fn atom_pair_default_adjacency_and_codes_match_smiles_topology() {
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
        let malformed_bond = BondEdge::new(1, 2, Bond::Single, None);
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

    #[test]
    fn rdk_fingerprint_graph_default_helper_builds_filtered_unique_bonds() {
        let malformed_bond = BondEdge::new(1, 2, Bond::Single, None);
        let graph = MalformedBondSmiles::new("CCO", 0, malformed_bond);

        let (atom_invariants, bonds, atom_bonds) = graph.rdk_fingerprint_graph(true, true);
        assert_eq!(atom_invariants, vec![10, 11, 12]);
        assert_eq!(
            bonds,
            vec![
                crate::traits::RdkFingerprintBond {
                    source: 0,
                    target: 1,
                    bond_code: 7,
                },
                crate::traits::RdkFingerprintBond {
                    source: 1,
                    target: 2,
                    bond_code: 7,
                },
            ]
        );
        assert_eq!(atom_bonds, vec![vec![0], vec![0, 1], vec![1]]);

        let (_atom_invariants, no_h_bonds, no_h_atom_bonds) =
            graph.rdk_fingerprint_graph(false, false);
        assert_eq!(
            no_h_bonds,
            vec![crate::traits::RdkFingerprintBond {
                source: 0,
                target: 1,
                bond_code: 1,
            }]
        );
        assert_eq!(no_h_atom_bonds, vec![vec![0], vec![0], vec![]]);
    }
}
