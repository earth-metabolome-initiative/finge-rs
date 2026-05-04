use core::hint::black_box;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use finge_rs::{
    CountEcfpFingerprint, EcfpFingerprint, Fingerprint, smiles_support::SmilesRdkitScratch,
};
use smiles_parser::smiles::Smiles;

mod common;

const BENCH_CASES: &[(u8, usize)] = &[
    (0, 64),
    (1, 128),
    (2, 128),
    (2, 2048),
    (4, 2048),
    (5, 2048),
    (5, 4096),
];

fn bench_corpus(c: &mut Criterion, corpus: &[Smiles]) {
    bench_bit_ecfp(c, corpus);
    bench_counted_ecfp(c, corpus);
}

fn bench_bit_ecfp(c: &mut Criterion, corpus: &[Smiles]) {
    let mut group = c.benchmark_group("ecfp_bits_smiles_with_rdkit_prep");
    common::configure_group(&mut group, corpus.len());

    for &(radius, fp_size) in BENCH_CASES {
        let fingerprint = EcfpFingerprint::new(radius, fp_size);
        group.bench_with_input(
            BenchmarkId::new(format!("r{radius}_n{fp_size}"), corpus.len()),
            &fingerprint,
            |b, fingerprint| {
                let mut scratch = SmilesRdkitScratch::default();
                b.iter(|| {
                    for graph in corpus {
                        let prepared = scratch.prepare(graph);
                        black_box(fingerprint.compute(&prepared));
                    }
                });
            },
        );
    }

    group.finish();
}

fn bench_counted_ecfp(c: &mut Criterion, corpus: &[Smiles]) {
    let mut group = c.benchmark_group("ecfp_counts_smiles_with_rdkit_prep");
    common::configure_group(&mut group, corpus.len());

    for &(radius, fp_size) in BENCH_CASES {
        let fingerprint = CountEcfpFingerprint::new(radius, fp_size);
        group.bench_with_input(
            BenchmarkId::new(format!("r{radius}_n{fp_size}"), corpus.len()),
            &fingerprint,
            |b, fingerprint| {
                let mut scratch = SmilesRdkitScratch::default();
                b.iter(|| {
                    for graph in corpus {
                        let prepared = scratch.prepare(graph);
                        black_box(fingerprint.compute(&prepared));
                    }
                });
            },
        );
    }

    group.finish();
}

fn ecfp_benchmarks(c: &mut Criterion) {
    let raw_corpus = common::load_smiles_corpus();
    bench_corpus(c, &raw_corpus);
}

criterion_group!(benches, ecfp_benchmarks);
criterion_main!(benches);
