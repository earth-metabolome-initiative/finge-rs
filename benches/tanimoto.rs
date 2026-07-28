use core::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use finge_rs::{
    BitFingerprint, EcfpFingerprint, Fingerprint, LshIndex, MinHash, SparseFingerprint,
    TanimotoIndex, TanimotoItem, smiles_support::SmilesRdkitScratch,
};
use smiles_parser::smiles::Smiles;

mod common;

/// ECFP radius and folded bit length for the benchmarked fingerprints.
const RADIUS: u8 = 2;
const FP_SIZE: usize = 2048;
/// Neighbours per query.
const QUERY_K: usize = 10;
/// MinHash permutations and band split for the LSH arm (rows = N / BANDS).
const N: usize = 512;
const BANDS: usize = 128;
/// Query sample size: `n << m`, a strided subsample of the `m`-molecule corpus.
const QUERY_SAMPLE: usize = 200;

/// Exact top-k by full pairwise Tanimoto, the no-index brute-force baseline.
/// Ordering matches the indices: descending similarity, ties ascending by id.
fn brute_topk<F: TanimotoItem>(corpus: &[F], query: &F, k: usize) -> Vec<(u64, f64)> {
    let mut scored: Vec<(u64, f64)> = corpus
        .iter()
        .enumerate()
        .map(|(id, item)| {
            (
                u64::try_from(id).expect("corpus fits u64"),
                query.tanimoto(item),
            )
        })
        .collect();
    scored.sort_by(|left, right| {
        right
            .1
            .partial_cmp(&left.1)
            .expect("scores are finite")
            .then(left.0.cmp(&right.0))
    });
    scored.truncate(k);
    scored
}

/// Folded-bit ECFP fingerprints plus their MinHash sketches (over the folded
/// active bins), one entry per corpus molecule.
fn dense_data(corpus: &[Smiles]) -> (Vec<BitFingerprint>, Vec<MinHash<u32, N>>) {
    let fingerprint = EcfpFingerprint::new(RADIUS, FP_SIZE);
    let mut scratch = SmilesRdkitScratch::default();
    let mut fingerprints = Vec::with_capacity(corpus.len());
    let mut sketches = Vec::with_capacity(corpus.len());
    for molecule in corpus {
        let graph = scratch.prepare(molecule);
        let bits = fingerprint.compute(&graph);
        let sketch: MinHash<u32, N> = bits
            .active_bits()
            .map(|bit| u64::try_from(bit).expect("bit index fits u64"))
            .collect();
        fingerprints.push(bits);
        sketches.push(sketch);
    }
    (fingerprints, sketches)
}

/// Unfolded ECFP fingerprints plus their MinHash sketches (over the raw Morgan
/// feature set), one entry per corpus molecule.
fn sparse_data(corpus: &[Smiles]) -> (Vec<SparseFingerprint<u32>>, Vec<MinHash<u32, N>>) {
    let fingerprint = EcfpFingerprint::new(RADIUS, FP_SIZE);
    let mut scratch = SmilesRdkitScratch::default();
    let mut fingerprints = Vec::with_capacity(corpus.len());
    let mut sketches = Vec::with_capacity(corpus.len());
    for molecule in corpus {
        let graph = scratch.prepare(molecule);
        fingerprints.push(fingerprint.unfolded(&graph));
        sketches.push(fingerprint.minhash::<_, u32, N>(&graph));
    }
    (fingerprints, sketches)
}

/// Compares the three retrieval methods for one representation: exact
/// brute-force (no index, zero build), banded MinHash LSH, and the exact
/// Tanimoto index. Build time and query time are measured in separate groups.
/// Fingerprints and sketches are precomputed, so timings exclude featurization.
fn bench_representation<F>(
    c: &mut Criterion,
    label: &str,
    fingerprints: &[F],
    sketches: &[MinHash<u32, N>],
) where
    F: TanimotoItem + Clone,
{
    let stride = (fingerprints.len() / QUERY_SAMPLE).max(1);
    let query_ids: Vec<usize> = (0..fingerprints.len()).step_by(stride).collect();

    // Build time. Brute force has no index, so its build time is zero and it has
    // no arm here.
    let mut build = c.benchmark_group(format!("{label}_build"));
    common::configure_group(&mut build, fingerprints.len());
    build.bench_function("lsh", |b| {
        b.iter(|| {
            let mut index = LshIndex::<u32, N, BANDS>::new();
            for sketch in sketches {
                index.insert(*sketch);
            }
            black_box(index)
        });
    });
    build.bench_function("tanimoto", |b| {
        b.iter(|| {
            let mut index = TanimotoIndex::new();
            for fingerprint in fingerprints {
                index.insert(fingerprint.clone());
            }
            index.freeze();
            black_box(index)
        });
    });
    build.finish();

    // Prebuild the indices once for the query group.
    let mut lsh = LshIndex::<u32, N, BANDS>::new();
    for sketch in sketches {
        lsh.insert(*sketch);
    }
    let mut tanimoto = TanimotoIndex::new();
    for fingerprint in fingerprints {
        tanimoto.insert(fingerprint.clone());
    }
    tanimoto.freeze();

    // Query time for the `n << m` sample against the full corpus.
    let mut query = c.benchmark_group(format!("{label}_query_k{QUERY_K}"));
    common::configure_group(&mut query, query_ids.len());
    query.bench_function("exact", |b| {
        b.iter(|| {
            for &id in &query_ids {
                black_box(brute_topk(fingerprints, &fingerprints[id], QUERY_K));
            }
        });
    });
    query.bench_function("lsh", |b| {
        b.iter(|| {
            for &id in &query_ids {
                black_box(lsh.query(&sketches[id], QUERY_K));
            }
        });
    });
    query.bench_function("tanimoto", |b| {
        b.iter(|| {
            for &id in &query_ids {
                black_box(tanimoto.query(&fingerprints[id], QUERY_K));
            }
        });
    });
    query.finish();
}

fn dense_benchmarks(c: &mut Criterion) {
    let corpus = common::load_smiles_corpus();
    let (fingerprints, sketches) = dense_data(&corpus);
    bench_representation(c, "dense", &fingerprints, &sketches);
}

fn sparse_benchmarks(c: &mut Criterion) {
    let corpus = common::load_smiles_corpus();
    let (fingerprints, sketches) = sparse_data(&corpus);
    bench_representation(c, "sparse", &fingerprints, &sketches);
}

criterion_group!(benches, dense_benchmarks, sparse_benchmarks);
criterion_main!(benches);
