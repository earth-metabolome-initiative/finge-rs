#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;

#[path = "common.rs"]
mod common;

#[derive(Debug, Arbitrary)]
struct SmilesInput {
    smiles: String,
}

fuzz_target!(|input: SmilesInput| {
    let Some(smiles) = common::parse_smiles(input.smiles) else {
        return;
    };
    common::fuzz_maccs(smiles);
});
