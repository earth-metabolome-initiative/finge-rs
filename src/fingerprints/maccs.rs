use alloc::boxed::Box;

use smarts_rs::{CompiledQuery, PreparedTarget};

use crate::{
    bit_fingerprint::BitFingerprint,
    fingerprint::Fingerprint,
    maccs_support::{
        MACCS_KEY_COUNT, MaccsBuildError, MaccsKeyDefinition, MaccsSpecialCase,
        compile_rdkit_maccs_queries, has_multiple_aromatic_rings, has_multiple_fragments,
        rdkit_maccs_keys,
    },
};

/// RDKit-style MACCS fingerprint built on top of `smarts-rs`.
#[derive(Debug, Clone)]
pub struct MaccsFingerprint {
    compiled_queries: Box<[Option<CompiledQuery>]>,
}

impl MaccsFingerprint {
    /// Builds the RDKit MACCS key matcher set.
    pub fn new() -> Result<Self, MaccsBuildError> {
        Ok(Self {
            compiled_queries: compile_rdkit_maccs_queries()?,
        })
    }

    fn key_is_set(&self, key: MaccsKeyDefinition, target: &PreparedTarget) -> bool {
        match key.special_case {
            Some(MaccsSpecialCase::Undefined) => false,
            Some(MaccsSpecialCase::AromaticRingCountGreaterThanOne) => {
                has_multiple_aromatic_rings(target)
            }
            Some(MaccsSpecialCase::FragmentCountGreaterThanOne) => has_multiple_fragments(target),
            None => {
                let compiled_query = self.compiled_queries[usize::from(key.key_id)]
                    .as_ref()
                    .expect("SMARTS-backed MACCS keys should compile");
                if key.count_threshold == 0 {
                    compiled_query.matches(target)
                } else {
                    compiled_query.match_count(target) > usize::from(key.count_threshold)
                }
            }
        }
    }
}

impl Fingerprint<PreparedTarget> for MaccsFingerprint {
    type Output = BitFingerprint;

    fn compute(&self, target: &PreparedTarget) -> Self::Output {
        let mut fingerprint = BitFingerprint::zeros(MACCS_KEY_COUNT);

        for key in rdkit_maccs_keys() {
            if self.key_is_set(key, target) {
                fingerprint.set(usize::from(key.key_id));
            }
        }

        fingerprint
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use smarts_rs::{CompiledQuery, PreparedTarget, QueryMol};
    use smiles_parser::{bond::Bond, smiles::Smiles};

    use super::MaccsFingerprint;
    use crate::{Fingerprint, test_fixtures::rdkit_maccs_fixture};

    const FIRST_RANDOM_1M_COUNTEREXAMPLE: &str = "CC[C@@H]1[C@@H](N1)C(=O)N(C)CC(=O)N(C)C(=C(C)C)C(=O)N[C@H]2CC3=CC(=CC(=C3)O)C4=CC5=C(C=C4)N(C(=C5CC(COC(=O)[C@@H]6CCCN(C2=O)N6)(C)C)C7=C(N=CC#C7)[C@H](C)OC)CC";

    fn prepare_target(smiles: &str) -> PreparedTarget {
        PreparedTarget::new(smiles.parse::<Smiles>().unwrap())
    }

    #[test]
    fn rdkit_maccs_known_examples_match() {
        let fingerprint = MaccsFingerprint::new().expect("MACCS SMARTS should compile");

        let cno = prepare_target("CNO");
        let ccc = prepare_target("CCC");

        assert_eq!(
            fingerprint.compute(&cno).active_bits().collect::<Vec<_>>(),
            alloc::vec![
                24, 68, 69, 71, 93, 94, 102, 124, 131, 139, 151, 158, 160, 161, 164
            ]
        );
        assert_eq!(
            fingerprint.compute(&ccc).active_bits().collect::<Vec<_>>(),
            alloc::vec![74, 114, 149, 155, 160]
        );
    }

    #[test]
    fn rdkit_maccs_fixture_matches_reference_corpus() {
        let fixture = rdkit_maccs_fixture();
        let fingerprint = MaccsFingerprint::new().expect("MACCS SMARTS should compile");

        for (smiles, expected_active_bits) in fixture.molecules.iter().zip(&fixture.active_bits) {
            let target = prepare_target(smiles);
            let observed_active_bits = fingerprint
                .compute(&target)
                .active_bits()
                .collect::<Vec<_>>();
            assert_eq!(
                observed_active_bits, *expected_active_bits,
                "MACCS mismatch for {smiles}"
            );
        }
    }

    #[test]
    fn maccs_key_17_is_set_for_simple_alkyne() {
        let target = prepare_target("CC#C");
        let fingerprint = MaccsFingerprint::new().expect("MACCS SMARTS should compile");
        let observed_active_bits = fingerprint
            .compute(&target)
            .active_bits()
            .collect::<Vec<_>>();

        assert!(observed_active_bits.contains(&17));
    }

    #[test]
    fn first_random_1m_counterexample_still_contains_a_carbon_carbon_triple_bond() {
        let smiles: Smiles = FIRST_RANDOM_1M_COUNTEREXAMPLE.parse().unwrap();

        let mut triple_cc_bonds = Vec::new();
        let mut atom_id = 0usize;
        while smiles.node_by_id(atom_id).is_some() {
            for edge in smiles.edges_for_node(atom_id) {
                let source = edge.source();
                let target = edge.target();
                if source >= target {
                    continue;
                }
                if edge.bond() == Bond::Triple
                    && smiles.node_by_id(source).and_then(|atom| atom.element())
                        == Some(elements_rs::Element::C)
                    && smiles.node_by_id(target).and_then(|atom| atom.element())
                        == Some(elements_rs::Element::C)
                {
                    triple_cc_bonds.push([source, target]);
                }
            }
            atom_id += 1;
        }

        assert_eq!(triple_cc_bonds, alloc::vec![[59, 60]]);
    }

    #[test]
    fn carbon_carbon_triple_bond_smarts_matches_first_random_1m_counterexample() {
        let target = prepare_target(FIRST_RANDOM_1M_COUNTEREXAMPLE);
        let query = "[#6]#[#6]".parse::<QueryMol>().unwrap();
        let compiled = CompiledQuery::new(query).unwrap();

        assert!(
            compiled.matches(&target),
            "prepared SMARTS matching should find the C#C bond in the counterexample"
        );
    }

    #[test]
    fn maccs_key_17_is_set_for_first_random_1m_counterexample() {
        let target = prepare_target(FIRST_RANDOM_1M_COUNTEREXAMPLE);
        let fingerprint = MaccsFingerprint::new().expect("MACCS SMARTS should compile");
        let observed_active_bits = fingerprint
            .compute(&target)
            .active_bits()
            .collect::<Vec<_>>();

        assert!(
            observed_active_bits.contains(&17),
            "MACCS key 17 ([#6]#[#6]) should be set for the counterexample, observed={observed_active_bits:?}"
        );
    }

    #[test]
    fn maccs_key_125_is_set_for_rdkit_aromatic_triple_bond_fused_ring() {
        let target = prepare_target("C1=CC2=CC#CC=C2C=C1");
        let fingerprint = MaccsFingerprint::new().expect("MACCS SMARTS should compile");
        let observed_active_bits = fingerprint
            .compute(&target)
            .active_bits()
            .collect::<Vec<_>>();

        assert!(
            observed_active_bits.contains(&125),
            "MACCS key 125 should be set for two RDKit-aromatic rings, observed={observed_active_bits:?}"
        );
    }
}
