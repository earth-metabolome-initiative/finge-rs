#![allow(dead_code)]

use finge_rs::{
    AtomPairFingerprint, BitFingerprint, CountEcfpFingerprint, CountFingerprint, EcfpFingerprint,
    Fingerprint, LayeredCountEcfpFingerprint, LayeredCountFingerprint, MaccsFingerprint,
    TopologicalTorsionFingerprint, smiles_support::SmilesRdkitScratch,
};
use smiles_parser::smiles::Smiles;
use smarts_validator::PreparedTarget;

const MAX_INPUT_BYTES: usize = 4096;

pub fn parse_smiles(input: String) -> Option<Smiles> {
    if input.is_empty() || input.len() > MAX_INPUT_BYTES {
        return None;
    }
    input.parse().ok()
}

pub fn prepare_graph<'a>(
    scratch: &'a mut SmilesRdkitScratch,
    smiles: &Smiles,
) -> Option<finge_rs::smiles_support::SmilesRdkitGraph<'a>> {
    scratch.try_prepare(smiles).ok()
}

pub fn assert_bit_fingerprint_basics(fingerprint: &BitFingerprint) {
    let active_bits = fingerprint.active_bits().collect::<Vec<_>>();
    assert!(active_bits.windows(2).all(|window| window[0] < window[1]));
    assert!(active_bits.iter().all(|&index| index < fingerprint.len()));
    assert!(active_bits.iter().all(|&index| fingerprint.contains(index)));
}

pub fn assert_count_fingerprint_basics(fingerprint: &CountFingerprint) {
    let active_counts = fingerprint.active_counts().collect::<Vec<_>>();
    assert!(active_counts.windows(2).all(|window| window[0].0 < window[1].0));
    assert!(active_counts
        .iter()
        .all(|&(index, count)| index < fingerprint.len() && count > 0));
    for &(index, count) in &active_counts {
        assert_eq!(fingerprint.count(index), count);
    }
}

pub fn assert_layered_count_fingerprint_basics(fingerprint: &LayeredCountFingerprint) {
    assert!(!fingerprint.is_empty());
    assert_eq!(fingerprint.formula(), &fingerprint.layers()[0]);
    for layer in fingerprint.layers() {
        assert_count_fingerprint_basics(layer);
    }
}

pub fn assert_count_matches_bit_presence(counts: &CountFingerprint, bits: &BitFingerprint) {
    let active_bits = bits.active_bits().collect::<Vec<_>>();
    let active_count_bits = counts.active_counts().map(|(index, _)| index).collect::<Vec<_>>();
    assert_eq!(active_bits, active_count_bits);
}

pub fn assert_layered_sums_match_total(layered: &LayeredCountFingerprint, total: &CountFingerprint) {
    assert_eq!(layered.formula().len(), total.len());
    for index in 0..total.len() {
        let layered_sum = layered
            .layers()
            .iter()
            .map(|layer| layer.count(index))
            .sum::<u32>();
        assert_eq!(layered_sum, total.count(index));
    }
}

pub fn fuzz_ecfp_family(smiles: Smiles) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };

    let bit = EcfpFingerprint::default().compute(&graph);
    let counted = CountEcfpFingerprint::default().compute(&graph);
    let layered = LayeredCountEcfpFingerprint::default().compute(&graph);
    let radius0 = CountEcfpFingerprint::new(0, 2048).compute(&graph);

    assert_bit_fingerprint_basics(&bit);
    assert_count_fingerprint_basics(&counted);
    assert_layered_count_fingerprint_basics(&layered);
    assert_count_matches_bit_presence(&counted, &bit);
    assert_layered_sums_match_total(&layered, &counted);
    assert_eq!(layered.formula(), &radius0);

    let bit_r1 = EcfpFingerprint::new(1, 2048).compute(&graph);
    let bit_r2 = EcfpFingerprint::new(2, 2048).compute(&graph);
    for index in bit_r1.active_bits() {
        assert!(bit_r2.contains(index));
    }

    let count_r1 = CountEcfpFingerprint::new(1, 2048).compute(&graph);
    let count_r2 = CountEcfpFingerprint::new(2, 2048).compute(&graph);
    for index in 0..count_r1.len() {
        assert!(count_r2.count(index) >= count_r1.count(index));
    }

    let explicit_h = smiles.with_explicit_hydrogens();
    let mut explicit_h_scratch = SmilesRdkitScratch::default();
    let Some(explicit_h_graph) = prepare_graph(&mut explicit_h_scratch, &explicit_h) else {
        return;
    };
    let explicit_h_bit = EcfpFingerprint::default().compute(&explicit_h_graph);
    let explicit_h_counted = CountEcfpFingerprint::default().compute(&explicit_h_graph);
    let explicit_h_layered = LayeredCountEcfpFingerprint::default().compute(&explicit_h_graph);

    assert_bit_fingerprint_basics(&explicit_h_bit);
    assert_count_fingerprint_basics(&explicit_h_counted);
    assert_layered_count_fingerprint_basics(&explicit_h_layered);
    assert_count_matches_bit_presence(&explicit_h_counted, &explicit_h_bit);
    assert_layered_sums_match_total(&explicit_h_layered, &explicit_h_counted);
}

pub fn fuzz_atom_pair(smiles: Smiles) {
    let raw_bits = AtomPairFingerprint::default().compute(&smiles);
    assert_bit_fingerprint_basics(&raw_bits);

    let without_count_sim = AtomPairFingerprint::default()
        .with_count_simulation(false)
        .compute(&smiles);
    assert_bit_fingerprint_basics(&without_count_sim);

    let explicit_h = smiles.with_explicit_hydrogens();
    let explicit_h_bits = AtomPairFingerprint::default().compute(&explicit_h);
    assert_bit_fingerprint_basics(&explicit_h_bits);

    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let prepared_bits = AtomPairFingerprint::default().compute(&graph);
    let prepared_without_count_sim = AtomPairFingerprint::default()
        .with_count_simulation(false)
        .compute(&graph);

    assert_bit_fingerprint_basics(&prepared_bits);
    assert_bit_fingerprint_basics(&prepared_without_count_sim);
}

pub fn fuzz_topological_torsion(smiles: Smiles) {
    let default_bits = TopologicalTorsionFingerprint::default().compute(&smiles);
    let no_count_sim = TopologicalTorsionFingerprint::default()
        .with_count_simulation(false)
        .compute(&smiles);
    let shortest_paths = TopologicalTorsionFingerprint::default()
        .with_only_shortest_paths(true)
        .compute(&smiles);
    let torsion5 = TopologicalTorsionFingerprint::default()
        .with_torsion_atom_count(5)
        .compute(&smiles);

    assert_bit_fingerprint_basics(&default_bits);
    assert_bit_fingerprint_basics(&no_count_sim);
    assert_bit_fingerprint_basics(&shortest_paths);
    assert_bit_fingerprint_basics(&torsion5);

    let explicit_h = smiles.with_explicit_hydrogens();
    let explicit_h_bits = TopologicalTorsionFingerprint::default().compute(&explicit_h);
    assert_bit_fingerprint_basics(&explicit_h_bits);

    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };

    let prepared_default = TopologicalTorsionFingerprint::default().compute(&graph);
    let prepared_no_count_sim = TopologicalTorsionFingerprint::default()
        .with_count_simulation(false)
        .compute(&graph);
    let prepared_shortest_paths = TopologicalTorsionFingerprint::default()
        .with_only_shortest_paths(true)
        .compute(&graph);
    let prepared_torsion5 = TopologicalTorsionFingerprint::default()
        .with_torsion_atom_count(5)
        .compute(&graph);

    assert_bit_fingerprint_basics(&prepared_default);
    assert_bit_fingerprint_basics(&prepared_no_count_sim);
    assert_bit_fingerprint_basics(&prepared_shortest_paths);
    assert_bit_fingerprint_basics(&prepared_torsion5);

    for index in prepared_shortest_paths.active_bits() {
        assert!(prepared_default.contains(index));
    }
}

pub fn fuzz_maccs(smiles: Smiles) {
    let fingerprint = MaccsFingerprint::new().expect("MACCS SMARTS should compile");
    let target = PreparedTarget::new(smiles.clone());
    let bits = fingerprint.compute(&target);
    assert_bit_fingerprint_basics(&bits);

    let explicit_h = smiles.with_explicit_hydrogens();
    let explicit_h_target = PreparedTarget::new(explicit_h);
    let explicit_h_bits = fingerprint.compute(&explicit_h_target);
    assert_bit_fingerprint_basics(&explicit_h_bits);

    let second_run = fingerprint.compute(&target);
    assert_eq!(
        bits.active_bits().collect::<Vec<_>>(),
        second_run.active_bits().collect::<Vec<_>>()
    );
}
