use core::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use finge_rs::{Fingerprint, MaccsFingerprint};
use smarts_rs::PreparedTarget;

mod common;

fn bench_corpus(c: &mut Criterion, corpus: &[PreparedTarget]) {
    let mut group = c.benchmark_group("maccs_smarts_prepared_target");
    common::configure_group(&mut group, corpus.len());

    let fingerprint = MaccsFingerprint::new().expect("MACCS SMARTS should compile");
    group.bench_function(format!("n{}", corpus.len()), |b| {
        b.iter(|| {
            for target in corpus {
                black_box(fingerprint.compute(target));
            }
        });
    });

    group.finish();
}

fn maccs_benchmarks(c: &mut Criterion) {
    let corpus = common::load_smiles_corpus()
        .into_iter()
        .map(PreparedTarget::new)
        .collect::<Vec<_>>();
    bench_corpus(c, &corpus);
}

criterion_group!(benches, maccs_benchmarks);
criterion_main!(benches);
