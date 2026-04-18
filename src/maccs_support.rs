use alloc::{boxed::Box, vec};
use core::fmt;

use smarts_parser::QueryMol;
use smarts_validator::{CompiledQuery, PreparedTarget, SmartsMatchError};

/// RDKit MACCS fingerprints use 167 bits with bit 0 left unused.
pub const MACCS_KEY_COUNT: usize = 167;

/// Keys that require counted substructure matches instead of plain boolean
/// presence.
pub const RDKIT_MACCS_THRESHOLD_KEY_IDS: [u16; 14] = [
    118, 120, 127, 130, 131, 136, 138, 140, 141, 142, 145, 146, 149, 159,
];

/// Non-SMARTS MACCS key cases from RDKit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaccsSpecialCase {
    /// Key 1 is intentionally undefined in RDKit.
    Undefined,
    /// Key 125 is set when the target contains more than one aromatic ring.
    AromaticRingCountGreaterThanOne,
    /// Key 166 is set when the target contains more than one connected
    /// component.
    FragmentCountGreaterThanOne,
}

/// One MACCS key definition taken directly from RDKit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaccsKeyDefinition {
    /// RDKit key id in the inclusive range `1..=166`.
    pub key_id: u16,
    /// SMARTS pattern when the key is SMARTS-backed.
    pub smarts: Option<&'static str>,
    /// Count threshold used by RDKit for `GetSubstructMatches`.
    pub count_threshold: u8,
    /// Non-SMARTS special-case behavior, when applicable.
    pub special_case: Option<MaccsSpecialCase>,
}

impl MaccsKeyDefinition {
    fn from_raw(key_id: u16, smarts: &'static str, count_threshold: u8) -> Self {
        let (smarts, special_case) = if smarts == "?" {
            let special_case = match key_id {
                1 => Some(MaccsSpecialCase::Undefined),
                125 => Some(MaccsSpecialCase::AromaticRingCountGreaterThanOne),
                166 => Some(MaccsSpecialCase::FragmentCountGreaterThanOne),
                _ => None,
            };
            (None, special_case)
        } else {
            (Some(smarts), None)
        };

        Self {
            key_id,
            smarts,
            count_threshold,
            special_case,
        }
    }
}

/// Build-time error while compiling the fixed RDKit MACCS SMARTS catalog.
#[derive(Debug)]
pub enum MaccsBuildError {
    /// One hard-coded SMARTS pattern failed to parse.
    Parse {
        key_id: u16,
        source: smarts_parser::SmartsParseError,
    },
    /// One parsed SMARTS pattern used validator features not yet supported by
    /// `smarts-validator`.
    Compile {
        key_id: u16,
        source: SmartsMatchError,
    },
}

impl fmt::Display for MaccsBuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse { key_id, source } => {
                write!(
                    f,
                    "failed to parse RDKit MACCS SMARTS for key {key_id}: {source}"
                )
            }
            Self::Compile { key_id, source } => {
                write!(
                    f,
                    "failed to compile RDKit MACCS SMARTS for key {key_id}: {source}"
                )
            }
        }
    }
}

/// Iterate the exact RDKit MACCS key catalog.
#[inline]
pub fn rdkit_maccs_keys() -> impl ExactSizeIterator<Item = MaccsKeyDefinition> + Clone {
    RDKIT_MACCS_KEY_DATA
        .iter()
        .copied()
        .map(|(key_id, smarts, count_threshold)| {
            MaccsKeyDefinition::from_raw(key_id, smarts, count_threshold)
        })
}

/// Compile every SMARTS-backed RDKit MACCS key through `smarts-rs`.
///
/// The returned slice is indexed directly by RDKit key id, so slot `0` is
/// intentionally empty.
pub fn compile_rdkit_maccs_queries() -> Result<Box<[Option<CompiledQuery>]>, MaccsBuildError> {
    let mut compiled = vec![None; MACCS_KEY_COUNT].into_boxed_slice();

    for key in rdkit_maccs_keys() {
        let Some(smarts) = key.smarts else {
            continue;
        };

        let query = smarts
            .parse::<QueryMol>()
            .map_err(|source| MaccsBuildError::Parse {
                key_id: key.key_id,
                source,
            })?;
        let compiled_query =
            CompiledQuery::new(query).map_err(|source| MaccsBuildError::Compile {
                key_id: key.key_id,
                source,
            })?;
        compiled[usize::from(key.key_id)] = Some(compiled_query);
    }

    Ok(compiled)
}

/// Returns whether the prepared target has more than one aromatic ring under
/// RDKit-default aromaticity.
pub fn has_multiple_aromatic_rings(target: &PreparedTarget) -> bool {
    let mut aromatic_ring_count = 0_u8;

    for cycle in target.target().symm_sssr_result().cycles() {
        if cycle.len() < 2 {
            continue;
        }

        let mut all_aromatic = true;
        for edge in cycle.windows(2) {
            if !is_aromatic_bond(target, edge[0], edge[1]) {
                all_aromatic = false;
                break;
            }
        }

        if all_aromatic && !is_aromatic_bond(target, cycle[cycle.len() - 1], cycle[0]) {
            all_aromatic = false;
        }

        if all_aromatic {
            aromatic_ring_count = aromatic_ring_count.saturating_add(1);
            if aromatic_ring_count > 1 {
                return true;
            }
        }
    }

    false
}

/// Returns whether the prepared target has more than one connected component.
pub fn has_multiple_fragments(target: &PreparedTarget) -> bool {
    let Some(first_component) = target.connected_component(0) else {
        return false;
    };

    (1..target.atom_count()).any(|atom_id| {
        target
            .connected_component(atom_id)
            .is_some_and(|component| component != first_component)
    })
}

fn is_aromatic_bond(target: &PreparedTarget, left_atom: usize, right_atom: usize) -> bool {
    target.bond(left_atom, right_atom) == Some(smarts_validator::BondLabel::Aromatic)
}

const RDKIT_MACCS_KEY_DATA: [(u16, &str, u8); 166] = [
    (1, "?", 0),
    (2, "[#104]", 0),
    (3, "[#32,#33,#34,#50,#51,#52,#82,#83,#84]", 0),
    (4, "[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]", 0),
    (5, "[Sc,Ti,Y,Zr,Hf]", 0),
    (6, "[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]", 0),
    (7, "[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]", 0),
    (8, "[!#6;!#1]1~*~*~*~1", 0),
    (9, "[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]", 0),
    (10, "[Be,Mg,Ca,Sr,Ba,Ra]", 0),
    (11, "*1~*~*~*~1", 0),
    (12, "[Cu,Zn,Ag,Cd,Au,Hg]", 0),
    (13, "[#8]~[#7](~[#6])~[#6]", 0),
    (14, "[#16]-[#16]", 0),
    (15, "[#8]~[#6](~[#8])~[#8]", 0),
    (16, "[!#6;!#1]1~*~*~1", 0),
    (17, "[#6]#[#6]", 0),
    (18, "[#5,#13,#31,#49,#81]", 0),
    (19, "*1~*~*~*~*~*~*~1", 0),
    (20, "[#14]", 0),
    (21, "[#6]=[#6](~[!#6;!#1])~[!#6;!#1]", 0),
    (22, "*1~*~*~1", 0),
    (23, "[#7]~[#6](~[#8])~[#8]", 0),
    (24, "[#7]-[#8]", 0),
    (25, "[#7]~[#6](~[#7])~[#7]", 0),
    (26, "[#6]=;@[#6](@*)@*", 0),
    (27, "[I]", 0),
    (28, "[!#6;!#1]~[CH2]~[!#6;!#1]", 0),
    (29, "[#15]", 0),
    (30, "[#6]~[!#6;!#1](~[#6])(~[#6])~*", 0),
    (31, "[!#6;!#1]~[F,Cl,Br,I]", 0),
    (32, "[#6]~[#16]~[#7]", 0),
    (33, "[#7]~[#16]", 0),
    (34, "[CH2]=*", 0),
    (35, "[Li,Na,K,Rb,Cs,Fr]", 0),
    (36, "[#16R]", 0),
    (37, "[#7]~[#6](~[#8])~[#7]", 0),
    (38, "[#7]~[#6](~[#6])~[#7]", 0),
    (39, "[#8]~[#16](~[#8])~[#8]", 0),
    (40, "[#16]-[#8]", 0),
    (41, "[#6]#[#7]", 0),
    (42, "F", 0),
    (43, "[!#6;!#1;!H0]~*~[!#6;!#1;!H0]", 0),
    (44, "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]", 0),
    (45, "[#6]=[#6]~[#7]", 0),
    (46, "Br", 0),
    (47, "[#16]~*~[#7]", 0),
    (48, "[#8]~[!#6;!#1](~[#8])(~[#8])", 0),
    (49, "[!+0]", 0),
    (50, "[#6]=[#6](~[#6])~[#6]", 0),
    (51, "[#6]~[#16]~[#8]", 0),
    (52, "[#7]~[#7]", 0),
    (53, "[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]", 0),
    (54, "[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]", 0),
    (55, "[#8]~[#16]~[#8]", 0),
    (56, "[#8]~[#7](~[#8])~[#6]", 0),
    (57, "[#8R]", 0),
    (58, "[!#6;!#1]~[#16]~[!#6;!#1]", 0),
    (59, "[#16]!:*:*", 0),
    (60, "[#16]=[#8]", 0),
    (61, "*~[#16](~*)~*", 0),
    (62, "*@*!@*@*", 0),
    (63, "[#7]=[#8]", 0),
    (64, "*@*!@[#16]", 0),
    (65, "c:n", 0),
    (66, "[#6]~[#6](~[#6])(~[#6])~*", 0),
    (67, "[!#6;!#1]~[#16]", 0),
    (68, "[!#6;!#1;!H0]~[!#6;!#1;!H0]", 0),
    (69, "[!#6;!#1]~[!#6;!#1;!H0]", 0),
    (70, "[!#6;!#1]~[#7]~[!#6;!#1]", 0),
    (71, "[#7]~[#8]", 0),
    (72, "[#8]~*~*~[#8]", 0),
    (73, "[#16]=*", 0),
    (74, "[CH3]~*~[CH3]", 0),
    (75, "*!@[#7]@*", 0),
    (76, "[#6]=[#6](~*)~*", 0),
    (77, "[#7]~*~[#7]", 0),
    (78, "[#6]=[#7]", 0),
    (79, "[#7]~*~*~[#7]", 0),
    (80, "[#7]~*~*~*~[#7]", 0),
    (81, "[#16]~*(~*)~*", 0),
    (82, "*~[CH2]~[!#6;!#1;!H0]", 0),
    (83, "[!#6;!#1]1~*~*~*~*~1", 0),
    (84, "[NH2]", 0),
    (85, "[#6]~[#7](~[#6])~[#6]", 0),
    (86, "[C;H2,H3][!#6;!#1][C;H2,H3]", 0),
    (87, "[F,Cl,Br,I]!@*@*", 0),
    (88, "[#16]", 0),
    (89, "[#8]~*~*~*~[#8]", 0),
    (
        90,
        "[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]",
        0,
    ),
    (
        91,
        "[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]",
        0,
    ),
    (92, "[#8]~[#6](~[#7])~[#6]", 0),
    (93, "[!#6;!#1]~[CH3]", 0),
    (94, "[!#6;!#1]~[#7]", 0),
    (95, "[#7]~*~*~[#8]", 0),
    (96, "*1~*~*~*~*~1", 0),
    (97, "[#7]~*~*~*~[#8]", 0),
    (98, "[!#6;!#1]1~*~*~*~*~*~1", 0),
    (99, "[#6]=[#6]", 0),
    (100, "*~[CH2]~[#7]", 0),
    (
        101,
        "[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]",
        0,
    ),
    (102, "[!#6;!#1]~[#8]", 0),
    (103, "Cl", 0),
    (104, "[!#6;!#1;!H0]~*~[CH2]~*", 0),
    (105, "*@*(@*)@*", 0),
    (106, "[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]", 0),
    (107, "[F,Cl,Br,I]~*(~*)~*", 0),
    (108, "[CH3]~*~*~*~[CH2]~*", 0),
    (109, "*~[CH2]~[#8]", 0),
    (110, "[#7]~[#6]~[#8]", 0),
    (111, "[#7]~*~[CH2]~*", 0),
    (112, "*~*(~*)(~*)~*", 0),
    (113, "[#8]!:*:*", 0),
    (114, "[CH3]~[CH2]~*", 0),
    (115, "[CH3]~*~[CH2]~*", 0),
    (116, "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]", 0),
    (117, "[#7]~*~[#8]", 0),
    (118, "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]", 1),
    (119, "[#7]=*", 0),
    (120, "[!#6;R]", 1),
    (121, "[#7;R]", 0),
    (122, "*~[#7](~*)~*", 0),
    (123, "[#8]~[#6]~[#8]", 0),
    (124, "[!#6;!#1]~[!#6;!#1]", 0),
    (125, "?", 0),
    (126, "*!@[#8]!@*", 0),
    (127, "*@*!@[#8]", 1),
    (
        128,
        "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]",
        0,
    ),
    (
        129,
        "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1)]",
        0,
    ),
    (130, "[!#6;!#1]~[!#6;!#1]", 1),
    (131, "[!#6;!#1;!H0]", 1),
    (132, "[#8]~*~[CH2]~*", 0),
    (133, "*@*!@[#7]", 0),
    (134, "[F,Cl,Br,I]", 0),
    (135, "[#7]!:*:*", 0),
    (136, "[#8]=*", 1),
    (137, "[!C;!c;R]", 0),
    (138, "[!#6;!#1]~[CH2]~*", 1),
    (139, "[O;!H0]", 0),
    (140, "[#8]", 3),
    (141, "[CH3]", 2),
    (142, "[#7]", 1),
    (143, "*@*!@[#8]", 0),
    (144, "*!:*:*!:*", 0),
    (145, "*1~*~*~*~*~*~1", 1),
    (146, "[#8]", 2),
    (147, "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]", 0),
    (148, "*~[!#6;!#1](~*)~*", 0),
    (149, "[C;H3,H4]", 1),
    (150, "*!@*@*!@*", 0),
    (151, "[#7;!H0]", 0),
    (152, "[#8]~[#6](~[#6])~[#6]", 0),
    (153, "[!#6;!#1]~[CH2]~*", 0),
    (154, "[#6]=[#8]", 0),
    (155, "*!@[CH2]!@*", 0),
    (156, "[#7]~*(~*)~*", 0),
    (157, "[#6]-[#8]", 0),
    (158, "[#6]-[#7]", 0),
    (159, "[#8]", 1),
    (160, "[C;H3,H4]", 0),
    (161, "[#7]", 0),
    (162, "a", 0),
    (163, "*1~*~*~*~*~*~1", 0),
    (164, "[#8]", 0),
    (165, "[R]", 0),
    (166, "?", 0),
];

#[cfg(test)]
mod tests {
    use smarts_validator::PreparedTarget;
    use smiles_parser::smiles::Smiles;

    use super::{
        MACCS_KEY_COUNT, MaccsKeyDefinition, MaccsSpecialCase, RDKIT_MACCS_THRESHOLD_KEY_IDS,
        compile_rdkit_maccs_queries, has_multiple_aromatic_rings, has_multiple_fragments,
        rdkit_maccs_keys,
    };
    use crate::test_fixtures::rdkit_maccs_fixture;

    #[test]
    fn rdkit_maccs_catalog_has_expected_shape() {
        let keys = rdkit_maccs_keys().collect::<alloc::vec::Vec<_>>();

        assert_eq!(keys.len(), MACCS_KEY_COUNT - 1);
        assert_eq!(
            keys[0],
            MaccsKeyDefinition {
                key_id: 1,
                smarts: None,
                count_threshold: 0,
                special_case: Some(MaccsSpecialCase::Undefined),
            }
        );
        assert_eq!(
            keys[124].special_case,
            Some(MaccsSpecialCase::AromaticRingCountGreaterThanOne)
        );
        assert_eq!(
            keys[165].special_case,
            Some(MaccsSpecialCase::FragmentCountGreaterThanOne)
        );
        assert_eq!(keys[113].smarts, Some("[CH3]~[CH2]~*"));
        assert_eq!(keys[117].count_threshold, 1);
    }

    #[test]
    fn rdkit_threshold_key_ids_match_source() {
        let threshold_key_ids = rdkit_maccs_keys()
            .filter(|key| key.count_threshold != 0)
            .map(|key| key.key_id)
            .collect::<alloc::vec::Vec<_>>();

        assert_eq!(threshold_key_ids, RDKIT_MACCS_THRESHOLD_KEY_IDS);
    }

    #[test]
    fn rdkit_maccs_queries_compile_through_smarts_rs() {
        let compiled = compile_rdkit_maccs_queries().expect("RDKit MACCS SMARTS should compile");

        assert_eq!(compiled.len(), MACCS_KEY_COUNT);
        assert!(compiled[0].is_none());
        assert!(compiled[1].is_none());
        assert!(compiled[2].is_some());
        assert!(compiled[125].is_none());
        assert!(compiled[166].is_none());

        let compiled_count = compiled.iter().filter(|entry| entry.is_some()).count();
        assert_eq!(compiled_count, MACCS_KEY_COUNT - 4);
    }

    #[test]
    fn aromatic_ring_special_case_matches_simple_examples() {
        let benzene = PreparedTarget::new("c1ccccc1".parse::<Smiles>().unwrap());
        let naphthalene = PreparedTarget::new("c1cccc2ccccc12".parse::<Smiles>().unwrap());

        assert!(!has_multiple_aromatic_rings(&benzene));
        assert!(has_multiple_aromatic_rings(&naphthalene));
    }

    #[test]
    fn fragment_special_case_matches_simple_examples() {
        let connected = PreparedTarget::new("CCO".parse::<Smiles>().unwrap());
        let disconnected = PreparedTarget::new("CC.O".parse::<Smiles>().unwrap());

        assert!(!has_multiple_fragments(&connected));
        assert!(has_multiple_fragments(&disconnected));
    }

    #[test]
    fn rdkit_maccs_fixture_loads_with_expected_shape() {
        let fixture = rdkit_maccs_fixture();

        assert_eq!(
            fixture.source.dataset,
            "scikit-fingerprints HIV test corpus"
        );
        assert_eq!(
            fixture.source.selection,
            "1024 parseable SMILES fixture in repo order"
        );
        assert_eq!(fixture.source.generator, "RDKit GetMACCSKeysFingerprint");
        assert_eq!(fixture.molecules.len(), 1_024);
        assert_eq!(fixture.active_bits.len(), fixture.molecules.len());
        assert!(
            fixture
                .active_bits
                .iter()
                .all(|active_bits| active_bits.iter().all(|&bit| bit < MACCS_KEY_COUNT))
        );
    }
}
