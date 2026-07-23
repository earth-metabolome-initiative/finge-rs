use core::hint::black_box;

use criterion::{
    BenchmarkGroup, BenchmarkId, Criterion, criterion_group, criterion_main, measurement::WallTime,
};
use finge_rs::{LshIndex, Map4Fingerprint, MinHash, smiles_support::SmilesRdkitScratch};

mod common;

/// MinHash permutation count for the benchmarked signatures.
const N: usize = 512;
/// Neighbours per query.
const QUERY_K: usize = 10;

fn signatures() -> Vec<MinHash<u32, N>> {
    let corpus = common::load_smiles_corpus();
    let fingerprint = Map4Fingerprint::default();
    let mut scratch = SmilesRdkitScratch::default();
    corpus
        .iter()
        .map(|molecule| {
            let graph = scratch.prepare(molecule);
            fingerprint.minhash::<_, u32, N>(&graph)
        })
        .collect()
}

fn build_index<const BANDS: usize>(signatures: &[MinHash<u32, N>]) -> LshIndex<u32, N, BANDS> {
    let mut index = LshIndex::<u32, N, BANDS>::new();
    for signature in signatures {
        index.insert(*signature);
    }
    index
}

// `BANDS` is a const generic now, so the old `for &bands in BANDS` sweep unrolls
// into one call per band count (`rows = N / BANDS`, the recall knob).
fn build_one<const BANDS: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    signatures: &[MinHash<u32, N>],
) {
    group.bench_function(BenchmarkId::new("bands", BANDS), |b| {
        b.iter(|| black_box(build_index::<BANDS>(signatures)));
    });
}

fn bench_build(c: &mut Criterion, signatures: &[MinHash<u32, N>]) {
    let mut group = c.benchmark_group("lsh_build");
    common::configure_group(&mut group, signatures.len());
    build_one::<128>(&mut group, signatures);
    build_one::<256>(&mut group, signatures);
    build_one::<512>(&mut group, signatures);
    group.finish();
}

fn query_one<const BANDS: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    signatures: &[MinHash<u32, N>],
    queries: &[MinHash<u32, N>],
) {
    let index = build_index::<BANDS>(signatures);
    group.bench_function(BenchmarkId::new("bands", BANDS), |b| {
        b.iter(|| {
            for query in queries {
                black_box(index.query(query, QUERY_K));
            }
        });
    });
}

fn bench_query(c: &mut Criterion, signatures: &[MinHash<u32, N>]) {
    let stride = (signatures.len() / 1000).max(1);
    let queries: Vec<MinHash<u32, N>> = signatures.iter().step_by(stride).copied().collect();
    let mut group = c.benchmark_group("lsh_query_k10");
    common::configure_group(&mut group, queries.len());
    query_one::<128>(&mut group, signatures, &queries);
    query_one::<256>(&mut group, signatures, &queries);
    query_one::<512>(&mut group, signatures, &queries);
    group.finish();
}

fn lsh_benchmarks(c: &mut Criterion) {
    let signatures = signatures();
    bench_build(c, &signatures);
    bench_query(c, &signatures);
}

criterion_group!(benches, lsh_benchmarks);
criterion_main!(benches);
