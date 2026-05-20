#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;

#[path = "common.rs"]
mod common;

#[derive(Debug, Arbitrary)]
struct MutatorInput {
    smiles: String,
    seed: u64,
}

fuzz_target!(|input: MutatorInput| {
    let Some(smiles) = common::parse_smiles(input.smiles) else {
        return;
    };
    common::fuzz_impossible_isotope_invariants(smiles, input.seed);
});
