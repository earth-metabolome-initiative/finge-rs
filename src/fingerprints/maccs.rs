use alloc::boxed::Box;

use smarts_validator::{PreparedTarget, match_count_compiled, matches_compiled};

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
    compiled_queries: Box<[Option<smarts_validator::CompiledQuery>]>,
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
                    matches_compiled(compiled_query, target)
                } else {
                    match_count_compiled(compiled_query, target) > usize::from(key.count_threshold)
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

    use smarts_validator::PreparedTarget;
    use smiles_parser::smiles::Smiles;

    use super::MaccsFingerprint;
    use crate::{Fingerprint, test_fixtures::rdkit_maccs_fixture};

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
}
