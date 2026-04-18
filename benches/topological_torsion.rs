use core::hint::black_box;
use std::io::Read;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use finge_rs::{Fingerprint, TopologicalTorsionFingerprint};
use flate2::read::GzDecoder;
use smiles_parser::smiles::Smiles;

const BENCH_SIZES: &[usize] = &[64, 128, 256, 512, 1024, 2048, 4096];

fn load_raw_corpus() -> Vec<Smiles> {
    let mut decoder = GzDecoder::new(
        &include_bytes!("../tests/fixtures/scikit_smallest_1024_parseable_smiles.txt.gz")[..],
    );
    let mut smiles = String::new();
    decoder
        .read_to_string(&mut smiles)
        .expect("benchmark fixture corpus should decompress");

    smiles
        .lines()
        .filter(|line| !line.is_empty())
        .map(str::to_owned)
        .map(|smiles| {
            smiles
                .parse()
                .expect("benchmark fixture SMILES should parse")
        })
        .collect()
}

fn bench_corpus(c: &mut Criterion, corpus: &[Smiles]) {
    let mut group = c.benchmark_group("topological_torsion_raw_smiles");
    group.sample_size(20);
    group.measurement_time(core::time::Duration::from_secs(3));
    group.warm_up_time(core::time::Duration::from_secs(1));
    group.throughput(Throughput::Elements(corpus.len() as u64));

    for &fp_size in BENCH_SIZES {
        let fingerprint = TopologicalTorsionFingerprint::new(fp_size);
        group.bench_with_input(
            BenchmarkId::new(format!("n{fp_size}"), corpus.len()),
            &fingerprint,
            |b, fingerprint| {
                b.iter(|| {
                    for graph in corpus {
                        black_box(fingerprint.compute(graph));
                    }
                });
            },
        );
    }

    group.finish();
}

fn topological_torsion_benchmarks(c: &mut Criterion) {
    let raw_corpus = load_raw_corpus();
    bench_corpus(c, &raw_corpus);
}

criterion_group!(benches, topological_torsion_benchmarks);
criterion_main!(benches);
