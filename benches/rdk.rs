use core::hint::black_box;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use finge_rs::{Fingerprint, RdkFingerprint, smiles_support::SmilesRdkitScratch};
use smiles_parser::smiles::Smiles;

mod common;

const BENCH_SIZES: &[usize] = &[64, 128, 256, 512, 1024, 2048, 4096];

fn bench_corpus(c: &mut Criterion, corpus: &[Smiles]) {
    let mut group = c.benchmark_group("rdk_smiles_with_rdkit_prep");
    common::configure_group(&mut group, corpus.len());

    for &fp_size in BENCH_SIZES {
        let fingerprint = RdkFingerprint::new(fp_size);
        group.bench_with_input(
            BenchmarkId::new(format!("n{fp_size}"), corpus.len()),
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

fn rdk_benchmarks(c: &mut Criterion) {
    let raw_corpus = common::load_smiles_corpus();
    bench_corpus(c, &raw_corpus);
}

criterion_group!(benches, rdk_benchmarks);
criterion_main!(benches);
