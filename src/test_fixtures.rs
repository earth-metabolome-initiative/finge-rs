use alloc::{string::String, vec::Vec};
use std::io::Read;
use std::sync::OnceLock;

use flate2::read::GzDecoder;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitEcfpSource {
    pub(crate) dataset: String,
    pub(crate) selection: String,
    pub(crate) generator: String,
    #[serde(rename = "includeChirality")]
    pub(crate) include_chirality: bool,
    #[serde(rename = "useBondTypes")]
    pub(crate) use_bond_types: bool,
    #[serde(rename = "includeRingMembership")]
    pub(crate) include_ring_membership: bool,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitEcfpCase {
    pub(crate) radius: u8,
    pub(crate) fp_size: usize,
    pub(crate) active_bits: Vec<Vec<usize>>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitEcfpFixture {
    pub(crate) source: RdkitEcfpSource,
    pub(crate) molecules: Vec<String>,
    pub(crate) cases: Vec<RdkitEcfpCase>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitCountedEcfpCase {
    pub(crate) radius: u8,
    pub(crate) fp_size: usize,
    pub(crate) active_counts: Vec<Vec<(usize, u32)>>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitCountedEcfpFixture {
    pub(crate) source: RdkitEcfpSource,
    pub(crate) molecules: Vec<String>,
    pub(crate) cases: Vec<RdkitCountedEcfpCase>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitLayeredCountedEcfpCase {
    pub(crate) radius: u8,
    pub(crate) fp_size: usize,
    pub(crate) layered_active_counts: Vec<Vec<Vec<(usize, u32)>>>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitLayeredCountedEcfpFixture {
    pub(crate) source: RdkitEcfpSource,
    pub(crate) molecules: Vec<String>,
    pub(crate) cases: Vec<RdkitLayeredCountedEcfpCase>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitAtomPairSource {
    pub(crate) dataset: String,
    pub(crate) selection: String,
    pub(crate) generator: String,
    #[serde(rename = "includeChirality")]
    pub(crate) include_chirality: bool,
    #[serde(rename = "countSimulation")]
    pub(crate) count_simulation: bool,
    #[serde(rename = "minDistance")]
    pub(crate) min_distance: u8,
    #[serde(rename = "maxDistance")]
    pub(crate) max_distance: u8,
    #[serde(rename = "use2D")]
    pub(crate) use_2d: bool,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitAtomPairCase {
    pub(crate) fp_size: usize,
    pub(crate) active_bits: Vec<Vec<usize>>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitAtomPairFixture {
    pub(crate) source: RdkitAtomPairSource,
    pub(crate) molecules: Vec<String>,
    pub(crate) cases: Vec<RdkitAtomPairCase>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitTopologicalTorsionSource {
    pub(crate) dataset: String,
    pub(crate) selection: String,
    pub(crate) generator: String,
    #[serde(rename = "includeChirality")]
    pub(crate) include_chirality: bool,
    #[serde(rename = "countSimulation")]
    pub(crate) count_simulation: bool,
    #[serde(rename = "onlyShortestPaths")]
    pub(crate) only_shortest_paths: bool,
    #[serde(rename = "torsionAtomCount")]
    pub(crate) torsion_atom_count: u8,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitTopologicalTorsionCase {
    pub(crate) fp_size: usize,
    pub(crate) active_bits: Vec<Vec<usize>>,
}

#[derive(Debug, Deserialize)]
pub(crate) struct RdkitTopologicalTorsionFixture {
    pub(crate) source: RdkitTopologicalTorsionSource,
    pub(crate) molecules: Vec<String>,
    pub(crate) cases: Vec<RdkitTopologicalTorsionCase>,
}

pub(crate) fn rdkit_ecfp_fixture() -> &'static RdkitEcfpFixture {
    static FIXTURE: OnceLock<RdkitEcfpFixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let mut decoder =
            GzDecoder::new(&include_bytes!("../tests/fixtures/rdkit_ecfp_reference.json.gz")[..]);
        let mut json = String::new();
        decoder
            .read_to_string(&mut json)
            .expect("RDKit ECFP reference fixture should decompress");
        serde_json::from_str(&json).expect("RDKit ECFP reference fixture should deserialize")
    })
}

pub(crate) fn rdkit_counted_ecfp_fixture() -> &'static RdkitCountedEcfpFixture {
    static FIXTURE: OnceLock<RdkitCountedEcfpFixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let mut decoder = GzDecoder::new(
            &include_bytes!("../tests/fixtures/rdkit_counted_ecfp_reference.json.gz")[..],
        );
        let mut json = String::new();
        decoder
            .read_to_string(&mut json)
            .expect("RDKit counted ECFP reference fixture should decompress");
        serde_json::from_str(&json)
            .expect("RDKit counted ECFP reference fixture should deserialize")
    })
}

pub(crate) fn rdkit_layered_counted_ecfp_fixture() -> &'static RdkitLayeredCountedEcfpFixture {
    static FIXTURE: OnceLock<RdkitLayeredCountedEcfpFixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let mut decoder = GzDecoder::new(
            &include_bytes!("../tests/fixtures/rdkit_layered_counted_ecfp_reference.json.gz")[..],
        );
        let mut json = String::new();
        decoder
            .read_to_string(&mut json)
            .expect("RDKit layered counted ECFP reference fixture should decompress");
        serde_json::from_str(&json)
            .expect("RDKit layered counted ECFP reference fixture should deserialize")
    })
}

pub(crate) fn rdkit_atom_pair_fixture() -> &'static RdkitAtomPairFixture {
    static FIXTURE: OnceLock<RdkitAtomPairFixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let mut decoder = GzDecoder::new(
            &include_bytes!("../tests/fixtures/rdkit_atom_pair_reference.json.gz")[..],
        );
        let mut json = String::new();
        decoder
            .read_to_string(&mut json)
            .expect("RDKit AtomPair reference fixture should decompress");
        serde_json::from_str(&json).expect("RDKit AtomPair reference fixture should deserialize")
    })
}

pub(crate) fn rdkit_topological_torsion_fixture() -> &'static RdkitTopologicalTorsionFixture {
    static FIXTURE: OnceLock<RdkitTopologicalTorsionFixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let mut decoder = GzDecoder::new(
            &include_bytes!("../tests/fixtures/rdkit_topological_torsion_reference.json.gz")[..],
        );
        let mut json = String::new();
        decoder
            .read_to_string(&mut json)
            .expect("RDKit Topological Torsion reference fixture should decompress");
        serde_json::from_str(&json)
            .expect("RDKit Topological Torsion reference fixture should deserialize")
    })
}
