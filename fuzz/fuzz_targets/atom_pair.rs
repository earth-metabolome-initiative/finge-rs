#![no_main]

mod common;

use libfuzzer_sys::fuzz_target;

fuzz_target!(|input: String| {
    let Some(smiles) = common::parse_smiles(input) else {
        return;
    };
    common::fuzz_atom_pair(smiles);
});
