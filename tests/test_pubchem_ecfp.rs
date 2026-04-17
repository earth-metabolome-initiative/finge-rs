#![cfg(feature = "smiles-support")]

use std::{
    env,
    fs::File,
    hint::black_box,
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
    time::{Duration, Instant},
};

use finge_rs::{EcfpFingerprint, Fingerprint, smiles_support::SmilesRdkitScratch};
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::{ThreadPoolBuilder, prelude::*};
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
    let threads = env_usize("PUBCHEM_ECFP_THREADS").unwrap_or_else(default_parallelism);
    let input_path = env::var_os("PUBCHEM_ECFP_INPUT").map(PathBuf::from);
    let fetch_options = DatasetFetchOptions {
        cache_dir: env::var_os("PUBCHEM_ECFP_CACHE_DIR").map(PathBuf::from),
        cache_mode: if env_flag("PUBCHEM_ECFP_REDOWNLOAD") {
            CacheMode::Redownload
        } else {
            CacheMode::UseCache
        },
        gzip_mode: GzipMode::Decompress,
    };
    let total_records = match input_path.as_deref() {
        Some(path) => limit.unwrap_or(count_smiles_records(path)?),
        None => limit.unwrap_or(PUBCHEM_CID_SMILES_RECORD_COUNT),
    };
    let progress_bar = pubchem_progress_bar(total_records);
    let started = Instant::now();
    let processed = AtomicUsize::new(0);
    let fingerprint = EcfpFingerprint::default();

    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(io::Error::other)?
        .install(|| match input_path.as_deref() {
            Some(path) => run_file_source(
                path,
                total_records,
                &fingerprint,
                &processed,
                &progress_bar,
                started,
            ),
            None => run_dataset_source(
                &fetch_options,
                total_records,
                &fingerprint,
                &processed,
                &progress_bar,
                started,
            ),
        })?;

    progress_bar.finish_with_message("ecfp complete");

    let processed = processed.load(Ordering::Relaxed);
    let elapsed = started.elapsed();
    eprintln!(
        "completed processed={processed} threads={threads} source={} elapsed={:.1}s rate={}/s",
        input_path
            .as_deref()
            .map_or_else(|| "pubchem".to_owned(), |path| path.display().to_string()),
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

fn run_dataset_source(
    fetch_options: &DatasetFetchOptions,
    total_records: usize,
    fingerprint: &EcfpFingerprint,
    processed: &AtomicUsize,
    progress_bar: &ProgressBar,
    started: Instant,
) -> io::Result<()> {
    PUBCHEM_SMILES
        .iter_smiles_with_options(fetch_options)
        .map_err(io::Error::other)?
        .take(total_records)
        .map(|line| line.map_err(|error| error.to_string()))
        .par_bridge()
        .try_for_each_init(SmilesRdkitScratch::default, |scratch, line| {
            process_smiles_line(scratch, fingerprint, processed, progress_bar, started, line)
        })
        .map_err(io::Error::other)
}

fn run_file_source(
    path: &Path,
    total_records: usize,
    fingerprint: &EcfpFingerprint,
    processed: &AtomicUsize,
    progress_bar: &ProgressBar,
    started: Instant,
) -> io::Result<()> {
    open_smiles_lines(path)?
        .take(total_records)
        .par_bridge()
        .try_for_each_init(SmilesRdkitScratch::default, |scratch, line| {
            process_smiles_line(scratch, fingerprint, processed, progress_bar, started, line)
        })
        .map_err(io::Error::other)
}

fn process_smiles_line(
    scratch: &mut SmilesRdkitScratch,
    fingerprint: &EcfpFingerprint,
    processed: &AtomicUsize,
    progress_bar: &ProgressBar,
    started: Instant,
    line: Result<String, String>,
) -> Result<(), String> {
    let smiles_text = line?;
    let smiles = smiles_text
        .parse::<Smiles>()
        .map_err(|error| format!("failed to parse SMILES:\n{}", error.render(&smiles_text)))?;

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
}

fn count_smiles_records(path: &Path) -> io::Result<usize> {
    open_smiles_lines(path)?.try_fold(0_usize, |count, line| {
        line.map(|_| count.saturating_add(1))
            .map_err(io::Error::other)
    })
}

fn open_smiles_lines(
    path: &Path,
) -> io::Result<Box<dyn Iterator<Item = Result<String, String>> + Send>> {
    let file = File::open(path)?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        let reader = BufReader::new(GzDecoder::new(file));
        Ok(Box::new(
            reader
                .lines()
                .map(|line| line.map_err(|error| error.to_string())),
        ))
    } else {
        let reader = BufReader::new(file);
        Ok(Box::new(
            reader
                .lines()
                .map(|line| line.map_err(|error| error.to_string())),
        ))
    }
}

fn default_parallelism() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
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
