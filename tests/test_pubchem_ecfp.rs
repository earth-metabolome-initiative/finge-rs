#![cfg(feature = "smiles-support")]

use std::{
    env,
    hint::black_box,
    io,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
    time::{Duration, Instant},
};

use finge_rs::{EcfpFingerprint, Fingerprint, smiles_support::SmilesEcfpScratch};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use smiles_parser::{
    datasets::{CacheMode, DatasetFetchOptions, GzipMode, PUBCHEM_SMILES, SmilesDatasetSource},
    smiles::Smiles,
};

const PUBCHEM_CID_SMILES_RECORD_COUNT: usize = 123_458_626;

#[test]
#[ignore = "This test downloads the full PubChem CID-SMILES corpus and fingerprints it in parallel."]
fn compute_ecfp_for_pubchem_smiles_corpus() -> Result<(), Box<dyn std::error::Error>> {
    // For a native-tuned run, use:
    // RUSTFLAGS='-C target-cpu=native' cargo test --release --features smiles-support --test test_pubchem_ecfp -- --ignored --nocapture
    require_release_build()?;

    let limit = env_usize("PUBCHEM_ECFP_LIMIT");
    let fetch_options = DatasetFetchOptions {
        cache_dir: env::var_os("PUBCHEM_ECFP_CACHE_DIR").map(PathBuf::from),
        cache_mode: if env_flag("PUBCHEM_ECFP_REDOWNLOAD") {
            CacheMode::Redownload
        } else {
            CacheMode::UseCache
        },
        gzip_mode: GzipMode::Decompress,
    };
    let total_records = limit.unwrap_or(PUBCHEM_CID_SMILES_RECORD_COUNT);
    let progress_bar = pubchem_progress_bar(total_records);
    let started = Instant::now();
    let processed = AtomicUsize::new(0);
    let fingerprint = EcfpFingerprint::default();

    PUBCHEM_SMILES
        .iter_smiles_with_options(&fetch_options)?
        .take(total_records)
        .par_bridge()
        .try_for_each_init(
            SmilesEcfpScratch::default,
            |scratch, line| -> Result<(), String> {
                let smiles_text = line.map_err(|error| error.to_string())?;
                let smiles = smiles_text.parse::<Smiles>().map_err(|error| {
                    format!("failed to parse SMILES:\n{}", error.render(&smiles_text))
                })?;

                let graph = scratch.prepare(&smiles);
                black_box(fingerprint.compute(&graph));

                let processed_now = processed.fetch_add(1, Ordering::Relaxed) + 1;
                progress_bar.inc(1);
                if processed_now % 100_000 == 0 {
                    progress_bar.set_message(format!(
                        "{}/s",
                        records_per_second(processed_now, started.elapsed())
                    ));
                }
                Ok(())
            },
        )
        .map_err(io::Error::other)?;

    progress_bar.finish_with_message("ecfp complete");

    let processed = processed.load(Ordering::Relaxed);
    let elapsed = started.elapsed();
    eprintln!(
        "completed processed={processed} elapsed={:.1}s rate={}/s",
        elapsed.as_secs_f64(),
        records_per_second(processed, elapsed),
    );

    Ok(())
}

fn pubchem_progress_bar(total_records: usize) -> ProgressBar {
    let progress_bar = ProgressBar::new(
        u64::try_from(total_records)
            .unwrap_or_else(|_| unreachable!("usize always fits into u64 on supported targets")),
    );
    progress_bar.enable_steady_tick(Duration::from_millis(100));
    progress_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({per_sec}, {eta}) {msg}",
        )
        .unwrap_or_else(|_| unreachable!("progress template is static and valid"))
        .progress_chars("=>-"),
    );
    progress_bar.set_message("ecfp");
    progress_bar
}

fn require_release_build() -> Result<(), io::Error> {
    if cfg!(debug_assertions) {
        return Err(io::Error::other(
            "run this ignored PubChem stress test with --release",
        ));
    }
    Ok(())
}

fn env_flag(name: &str) -> bool {
    env::var_os(name).is_some_and(|value| value != "0")
}

fn env_usize(name: &str) -> Option<usize> {
    env::var(name).ok().and_then(|raw| raw.parse().ok())
}

fn records_per_second(records: usize, elapsed: Duration) -> usize {
    let nanos = elapsed.as_nanos();
    if nanos == 0 {
        return 0;
    }

    let records =
        u128::try_from(records).unwrap_or_else(|_| unreachable!("usize always fits into u128"));
    let per_second = (records
        .saturating_mul(1_000_000_000)
        .saturating_add(nanos / 2))
        / nanos;
    usize::try_from(per_second).unwrap_or(usize::MAX)
}
