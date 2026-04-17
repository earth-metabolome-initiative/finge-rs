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
