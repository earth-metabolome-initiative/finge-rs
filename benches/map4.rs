use core::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use finge_rs::Map4Fingerprint;

mod common;

/// MinHash permutation count for the sketched signatures.
const N: usize = 512;

/// Times MAP4 MinHash sketch construction over the whole corpus. This is the
/// shingle-to-sketch path `benches/lsh.rs` builds in its untimed setup,
/// isolated here so the shingle enumeration and MinHash insertion cost is
/// measured on its own. Sketching runs over the parsed `Smiles` directly, which
/// yields the same shingle set as the RDKit-normalized graph
/// (`map4_normalized_path_matches_raw`).
fn bench_sketch(c: &mut Criterion) {
    let corpus = common::load_smiles_corpus();
    let fingerprint = Map4Fingerprint::default();

    let mut group = c.benchmark_group("map4_sketch");
    common::configure_group(&mut group, corpus.len());
    group.bench_function("minhash_u32_512", |b| {
        b.iter(|| {
            for molecule in &corpus {
                black_box(fingerprint.minhash::<_, u32, N>(molecule));
            }
        });
    });
    group.finish();
}

criterion_group!(benches, bench_sketch);
criterion_main!(benches);
