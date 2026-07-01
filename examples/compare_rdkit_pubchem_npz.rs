use std::{
    collections::{BTreeMap, BTreeSet},
    env,
    error::Error,
    fmt,
    fs::{self, File},
    io::{BufWriter, Cursor, Read, Write},
    path::{Path, PathBuf},
};

#[cfg(feature = "smarts-support")]
use finge_rs::MaccsFingerprint;
use finge_rs::{
    AtomPairFingerprint, CountEcfpFingerprint, EcfpFingerprint, Fingerprint,
    LayeredCountEcfpFingerprint, RdkFingerprint, SmilesRdkitScratch, TopologicalTorsionFingerprint,
};
use rayon::prelude::*;
#[cfg(feature = "smarts-support")]
use smarts_rs::PreparedTarget;
use smiles_parser::smiles::Smiles;
use zip::ZipArchive;

type DynError = Box<dyn Error + Send + Sync>;

#[cfg(feature = "smarts-support")]
type MaccsContext = MaccsFingerprint;

#[cfg(not(feature = "smarts-support"))]
#[derive(Debug, Default)]
struct MaccsContext;

#[derive(Debug)]
struct Args {
    input: PathBuf,
    limit: Option<usize>,
    max_shards: Option<usize>,
    mismatch_limit: usize,
    mismatch_output: Option<PathBuf>,
    threads: Option<usize>,
    fingerprints: BTreeSet<String>,
    skip_unparseable: bool,
}

#[derive(Debug)]
struct Shard {
    path: PathBuf,
    smiles: Vec<String>,
    bit_arrays: BTreeMap<String, U8Array>,
    count_arrays: BTreeMap<String, CsrArray>,
}

#[derive(Debug)]
struct U8Array {
    shape: Vec<usize>,
    data: Vec<u8>,
}

#[derive(Debug)]
struct CsrArray {
    indptr: Vec<u64>,
    indices: Vec<u32>,
    counts: Vec<u32>,
}

#[derive(Debug, Clone)]
enum BitCase {
    Ecfp {
        key: String,
        radius: u8,
        fp_size: usize,
    },
    AtomPair {
        key: String,
        fp_size: usize,
    },
    TopologicalTorsion {
        key: String,
        fp_size: usize,
    },
    Maccs {
        key: String,
    },
    Rdk {
        key: String,
        fp_size: usize,
    },
}

#[derive(Debug, Clone)]
enum CountCase {
    Ecfp {
        key: String,
        radius: u8,
        fp_size: usize,
    },
    LayeredEcfp {
        key: String,
        radius: u8,
        layer: usize,
        fp_size: usize,
    },
}

#[derive(Debug, Default)]
struct CaseSummary {
    checked: usize,
    mismatches: usize,
    samples: Vec<MismatchSample>,
}

#[derive(Debug)]
struct MismatchSample {
    row: usize,
    smiles: String,
    expected: String,
    observed: String,
}

#[derive(Debug)]
struct RowMismatch {
    key: String,
    row: usize,
    smiles: String,
    expected: String,
    observed: String,
}

#[derive(Debug)]
enum RowComparison {
    Compared(Vec<RowMismatch>),
    SkippedParse {
        row: usize,
        smiles: String,
        error: String,
    },
}

#[derive(Debug)]
struct SkippedSample {
    row: usize,
    smiles: String,
    error: String,
}

struct ComparisonContext<'a> {
    bit_cases: &'a [BitCase],
    count_cases: &'a [CountCase],
    maccs: &'a MaccsContext,
    skip_unparseable: bool,
}

impl Args {
    fn parse() -> Result<Self, DynError> {
        let mut args = env::args().skip(1);
        let mut input = None;
        let mut limit = None;
        let mut max_shards = None;
        let mut mismatch_limit = 8usize;
        let mut mismatch_output = None;
        let mut threads = None;
        let mut fingerprints = BTreeSet::new();
        let mut skip_unparseable = false;

        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--input-dir" | "--input" => {
                    input = args.next().map(PathBuf::from);
                }
                "--limit" => {
                    limit = Some(next_arg(&mut args, "--limit")?.parse()?);
                }
                "--max-shards" => {
                    max_shards = Some(next_arg(&mut args, "--max-shards")?.parse()?);
                }
                "--mismatch-limit" => {
                    mismatch_limit = next_arg(&mut args, "--mismatch-limit")?.parse()?;
                }
                "--mismatch-output" => {
                    mismatch_output =
                        Some(PathBuf::from(next_arg(&mut args, "--mismatch-output")?));
                }
                "--threads" => {
                    threads = Some(next_arg(&mut args, "--threads")?.parse()?);
                }
                "--fingerprint" => {
                    fingerprints.insert(next_arg(&mut args, "--fingerprint")?);
                }
                "--skip-unparseable" => {
                    skip_unparseable = true;
                }
                "--help" | "-h" => {
                    print_usage();
                    std::process::exit(0);
                }
                other => return Err(format!("unknown argument: {other}").into()),
            }
        }

        Ok(Self {
            input: input.ok_or("--input-dir is required")?,
            limit,
            max_shards,
            mismatch_limit,
            mismatch_output,
            threads,
            fingerprints,
            skip_unparseable,
        })
    }

    fn wants_family(&self, family: &str) -> bool {
        self.fingerprints.is_empty() || self.fingerprints.contains(family)
    }
}

fn next_arg(args: &mut impl Iterator<Item = String>, name: &str) -> Result<String, DynError> {
    args.next()
        .ok_or_else(|| format!("{name} requires a value").into())
}

fn print_usage() {
    eprintln!(
        "usage: cargo run --release --example compare_rdkit_pubchem_npz -- --input-dir /tmp/finge_rs_rdkit_pubchem_default [--threads 64] [--limit N] [--max-shards N] [--skip-unparseable] [--mismatch-output failures.tsv] [--fingerprint ecfp|counted-ecfp|layered-counted-ecfp|atom-pair|topological-torsion|maccs|rdk]"
    );
}

fn main() -> Result<(), DynError> {
    let args = Args::parse()?;
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    let mut paths = npz_paths(&args.input)?;
    if let Some(max_shards) = args.max_shards {
        paths.truncate(max_shards);
    }
    if paths.is_empty() {
        return Err(format!("no .npz shards found at {}", args.input.display()).into());
    }

    let mut summaries = BTreeMap::<String, CaseSummary>::new();
    let mut global_row = 0usize;
    let mut compared_rows = 0usize;
    let mut skipped_rows = 0usize;
    let mut skipped_samples = Vec::new();
    let mut remaining = args.limit;
    let mut mismatch_writer = args
        .mismatch_output
        .as_ref()
        .map(|path| -> Result<BufWriter<File>, DynError> {
            let mut writer = BufWriter::new(File::create(path)?);
            writeln!(writer, "key\trow\tsmiles\texpected\tobserved")?;
            Ok(writer)
        })
        .transpose()?;
    let maccs = build_maccs_context(&args)?;

    for (shard_index, path) in paths.iter().enumerate() {
        if remaining == Some(0) {
            break;
        }

        eprintln!("reading shard {}: {}", shard_index + 1, path.display());
        let shard = Shard::open(path, &args)?;
        let bit_cases = shard.bit_cases(&args);
        let count_cases = shard.count_cases(&args);
        let row_limit = remaining.map_or(shard.smiles.len(), |limit| limit.min(shard.smiles.len()));
        let context = ComparisonContext {
            bit_cases: &bit_cases,
            count_cases: &count_cases,
            maccs: &maccs,
            skip_unparseable: args.skip_unparseable,
        };

        let row_mismatches = (0..row_limit)
            .into_par_iter()
            .map_init(SmilesRdkitScratch::default, |scratch, row| {
                compare_row(&shard, row, global_row + row, &context, scratch)
            })
            .collect::<Vec<_>>();

        for row_result in row_mismatches {
            match row_result? {
                RowComparison::Compared(mismatches) => {
                    compared_rows += 1;
                    for case in &bit_cases {
                        summaries.entry(case.key().to_owned()).or_default().checked += 1;
                    }
                    for case in &count_cases {
                        summaries.entry(case.key().to_owned()).or_default().checked += 1;
                    }
                    for mismatch in mismatches {
                        if let Some(writer) = mismatch_writer.as_mut() {
                            write_mismatch_tsv(writer, &mismatch)?;
                        }
                        summaries
                            .entry(mismatch.key.clone())
                            .or_default()
                            .record_formatted(mismatch, args.mismatch_limit);
                    }
                }
                RowComparison::SkippedParse { row, smiles, error } => {
                    skipped_rows += 1;
                    if skipped_samples.len() < args.mismatch_limit {
                        skipped_samples.push(SkippedSample { row, smiles, error });
                    }
                }
            }
        }
        global_row += row_limit;

        if let Some(remaining_rows) = remaining.as_mut() {
            *remaining_rows -= row_limit;
        }
    }

    if let Some(writer) = mismatch_writer.as_mut() {
        writer.flush()?;
    }
    print_summary(&summaries);
    print_skipped_summary(skipped_rows, &skipped_samples);
    let total_mismatches = summaries
        .values()
        .map(|summary| summary.mismatches)
        .sum::<usize>();
    if total_mismatches == 0 {
        println!("all selected fingerprints matched RDKit for {compared_rows} parseable molecules");
        Ok(())
    } else {
        Err(format!("found {total_mismatches} total fingerprint mismatches").into())
    }
}

#[cfg(feature = "smarts-support")]
fn build_maccs_context(_args: &Args) -> Result<MaccsContext, DynError> {
    MaccsFingerprint::new()
        .map_err(|error| format!("failed to build MACCS fingerprint: {error}").into())
}

#[cfg(not(feature = "smarts-support"))]
fn build_maccs_context(args: &Args) -> Result<MaccsContext, DynError> {
    if args.wants_family("maccs") {
        Err("MACCS comparison requires the smarts-support feature".into())
    } else {
        Ok(MaccsContext)
    }
}

fn write_mismatch_tsv(writer: &mut impl Write, mismatch: &RowMismatch) -> Result<(), DynError> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}",
        escape_tsv(&mismatch.key),
        mismatch.row,
        escape_tsv(&mismatch.smiles),
        escape_tsv(&mismatch.expected),
        escape_tsv(&mismatch.observed),
    )?;
    Ok(())
}

fn escape_tsv(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
}

fn npz_paths(input: &Path) -> Result<Vec<PathBuf>, DynError> {
    if input.is_file() {
        return Ok(vec![input.to_owned()]);
    }

    let mut paths = fs::read_dir(input)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<Result<Vec<_>, _>>()?;
    paths.retain(|path| path.extension().is_some_and(|extension| extension == "npz"));
    paths.sort();
    Ok(paths)
}

impl Shard {
    fn open(path: &Path, args: &Args) -> Result<Self, DynError> {
        let file = File::open(path)?;
        let mut archive = ZipArchive::new(file)?;
        let names = archive_names(&mut archive)?;

        let smiles_offsets = read_npy_u64(&mut archive, "smiles_offsets.npy")?;
        let smiles_utf8 = read_npy_u8(&mut archive, "smiles_utf8.npy")?;
        let smiles = decode_smiles(&smiles_offsets, &smiles_utf8.data)?;

        let mut bit_arrays = BTreeMap::new();
        let mut count_keys = BTreeSet::new();
        for name in &names {
            if let Some(key) = name.strip_suffix("__bits.npy") {
                if parse_bit_case(key).is_some_and(|case| args.wants_family(case.family())) {
                    bit_arrays.insert(key.to_owned(), read_npy_u8(&mut archive, name)?);
                }
            } else if let Some(key) = name.strip_suffix("__indptr.npy") {
                if parse_count_case(key).is_some_and(|case| args.wants_family(case.family())) {
                    count_keys.insert(key.to_owned());
                }
            }
        }

        let mut count_arrays = BTreeMap::new();
        for key in count_keys {
            count_arrays.insert(
                key.clone(),
                CsrArray {
                    indptr: read_npy_u64(&mut archive, &format!("{key}__indptr.npy"))?,
                    indices: read_npy_u32(&mut archive, &format!("{key}__indices.npy"))?,
                    counts: read_npy_u32(&mut archive, &format!("{key}__counts.npy"))?,
                },
            );
        }

        Ok(Self {
            path: path.to_owned(),
            smiles,
            bit_arrays,
            count_arrays,
        })
    }

    fn bit_cases(&self, args: &Args) -> Vec<BitCase> {
        self.bit_arrays
            .keys()
            .filter_map(|key| parse_bit_case(key))
            .filter(|case| args.wants_family(case.family()))
            .collect()
    }

    fn count_cases(&self, args: &Args) -> Vec<CountCase> {
        self.count_arrays
            .keys()
            .filter_map(|key| parse_count_case(key))
            .filter(|case| args.wants_family(case.family()))
            .collect()
    }

    fn bit_row(&self, key: &str, row: usize) -> Result<&[u8], DynError> {
        let array = self
            .bit_arrays
            .get(key)
            .ok_or_else(|| format!("missing bit array {key} in {}", self.path.display()))?;
        if array.shape.len() != 2 {
            return Err(format!("bit array {key} should be rank-2").into());
        }
        let row_count = array.shape[0];
        let row_width = array.shape[1];
        if row >= row_count {
            return Err(format!("row {row} out of bounds for bit array {key}").into());
        }
        let start = row * row_width;
        Ok(&array.data[start..start + row_width])
    }

    fn count_row(&self, key: &str, row: usize) -> Result<Vec<(usize, u32)>, DynError> {
        let array = self
            .count_arrays
            .get(key)
            .ok_or_else(|| format!("missing count array {key} in {}", self.path.display()))?;
        let start = usize::try_from(array.indptr[row])?;
        let end = usize::try_from(array.indptr[row + 1])?;
        Ok((start..end)
            .map(|index| (array.indices[index] as usize, array.counts[index]))
            .collect())
    }
}

fn archive_names<R: Read + std::io::Seek>(
    archive: &mut ZipArchive<R>,
) -> Result<Vec<String>, DynError> {
    let mut names = Vec::with_capacity(archive.len());
    for index in 0..archive.len() {
        names.push(archive.by_index(index)?.name().to_owned());
    }
    Ok(names)
}

impl BitCase {
    fn key(&self) -> &str {
        match self {
            Self::Ecfp { key, .. }
            | Self::AtomPair { key, .. }
            | Self::TopologicalTorsion { key, .. }
            | Self::Maccs { key }
            | Self::Rdk { key, .. } => key,
        }
    }

    fn family(&self) -> &str {
        match self {
            Self::Ecfp { .. } => "ecfp",
            Self::AtomPair { .. } => "atom-pair",
            Self::TopologicalTorsion { .. } => "topological-torsion",
            Self::Maccs { .. } => "maccs",
            Self::Rdk { .. } => "rdk",
        }
    }

    fn fp_size(&self) -> usize {
        match self {
            Self::Ecfp { fp_size, .. }
            | Self::AtomPair { fp_size, .. }
            | Self::TopologicalTorsion { fp_size, .. }
            | Self::Rdk { fp_size, .. } => *fp_size,
            Self::Maccs { .. } => 167,
        }
    }

    fn compute(
        &self,
        smiles: &Smiles,
        graph: &finge_rs::SmilesRdkitGraph<'_>,
        maccs: &MaccsContext,
    ) -> Result<Vec<u8>, DynError> {
        let bit_count = self.fp_size();
        let mut packed = vec![0_u8; byte_width(bit_count)];
        match self {
            Self::Ecfp {
                radius, fp_size, ..
            } => {
                set_packed_bits(
                    &mut packed,
                    EcfpFingerprint::new(*radius, *fp_size)
                        .compute(graph)
                        .active_bits(),
                );
            }
            Self::AtomPair { fp_size, .. } => {
                set_packed_bits(
                    &mut packed,
                    AtomPairFingerprint::new(*fp_size)
                        .compute(graph)
                        .active_bits(),
                );
            }
            Self::TopologicalTorsion { fp_size, .. } => {
                set_packed_bits(
                    &mut packed,
                    TopologicalTorsionFingerprint::new(*fp_size)
                        .compute(graph)
                        .active_bits(),
                );
            }
            Self::Maccs { .. } => {
                #[cfg(feature = "smarts-support")]
                {
                    let target = PreparedTarget::new(smiles.clone());
                    set_packed_bits(&mut packed, maccs.compute(&target).active_bits());
                }
                #[cfg(not(feature = "smarts-support"))]
                {
                    let _ = (smiles, maccs);
                    return Err("MACCS comparison requires the smarts-support feature".into());
                }
            }
            Self::Rdk { fp_size, .. } => {
                set_packed_bits(
                    &mut packed,
                    RdkFingerprint::new(*fp_size).compute(graph).active_bits(),
                );
            }
        }
        Ok(packed)
    }
}

impl CountCase {
    fn key(&self) -> &str {
        match self {
            Self::Ecfp { key, .. } | Self::LayeredEcfp { key, .. } => key,
        }
    }

    fn family(&self) -> &str {
        match self {
            Self::Ecfp { .. } => "counted-ecfp",
            Self::LayeredEcfp { .. } => "layered-counted-ecfp",
        }
    }

    fn compute(&self, graph: &finge_rs::SmilesRdkitGraph<'_>) -> Vec<(usize, u32)> {
        match self {
            Self::Ecfp {
                radius, fp_size, ..
            } => CountEcfpFingerprint::new(*radius, *fp_size)
                .compute(graph)
                .active_counts()
                .collect(),
            Self::LayeredEcfp {
                radius,
                layer,
                fp_size,
                ..
            } => LayeredCountEcfpFingerprint::new(*radius, *fp_size)
                .compute(graph)
                .layer(*layer)
                .expect("generated layer case should be in range")
                .active_counts()
                .collect(),
        }
    }
}

fn parse_bit_case(key: &str) -> Option<BitCase> {
    if let Some((radius, fp_size)) = parse_radius_size(key, "ecfp_bits_r") {
        return Some(BitCase::Ecfp {
            key: key.to_owned(),
            radius,
            fp_size,
        });
    }
    if let Some(fp_size) = parse_size(key, "atom_pair_bits_n") {
        return Some(BitCase::AtomPair {
            key: key.to_owned(),
            fp_size,
        });
    }
    if let Some(fp_size) = parse_size(key, "topological_torsion_bits_n") {
        return Some(BitCase::TopologicalTorsion {
            key: key.to_owned(),
            fp_size,
        });
    }
    if key == "maccs_bits" {
        return Some(BitCase::Maccs {
            key: key.to_owned(),
        });
    }
    if let Some(fp_size) = parse_size(key, "rdk_bits_n") {
        return Some(BitCase::Rdk {
            key: key.to_owned(),
            fp_size,
        });
    }
    None
}

fn parse_count_case(key: &str) -> Option<CountCase> {
    if let Some((radius, fp_size)) = parse_radius_size(key, "ecfp_counts_r") {
        return Some(CountCase::Ecfp {
            key: key.to_owned(),
            radius,
            fp_size,
        });
    }

    let rest = key.strip_prefix("ecfp_layer_counts_r")?;
    let (radius, rest) = rest.split_once("_layer")?;
    let (layer, fp_size) = rest.split_once("_n")?;
    Some(CountCase::LayeredEcfp {
        key: key.to_owned(),
        radius: radius.parse().ok()?,
        layer: layer.parse().ok()?,
        fp_size: fp_size.parse().ok()?,
    })
}

fn parse_radius_size(key: &str, prefix: &str) -> Option<(u8, usize)> {
    let rest = key.strip_prefix(prefix)?;
    let (radius, fp_size) = rest.split_once("_n")?;
    Some((radius.parse().ok()?, fp_size.parse().ok()?))
}

fn parse_size(key: &str, prefix: &str) -> Option<usize> {
    key.strip_prefix(prefix)?.parse().ok()
}

impl CaseSummary {
    fn record_formatted(&mut self, mismatch: RowMismatch, mismatch_limit: usize) {
        self.mismatches += 1;
        if self.samples.len() < mismatch_limit {
            self.samples.push(MismatchSample {
                row: mismatch.row,
                smiles: mismatch.smiles,
                expected: mismatch.expected,
                observed: mismatch.observed,
            });
        }
    }
}

fn compare_row(
    shard: &Shard,
    row: usize,
    global_row: usize,
    context: &ComparisonContext<'_>,
    scratch: &mut SmilesRdkitScratch,
) -> Result<RowComparison, DynError> {
    let smiles_text = &shard.smiles[row];
    let smiles: Smiles = match smiles_text.parse() {
        Ok(smiles) => smiles,
        Err(error) if context.skip_unparseable => {
            return Ok(RowComparison::SkippedParse {
                row: global_row,
                smiles: smiles_text.clone(),
                error: error.to_string(),
            });
        }
        Err(error) => {
            return Err(format!("failed to parse row {global_row}: {smiles_text}: {error}").into());
        }
    };
    let graph = match scratch.try_prepare(&smiles) {
        Ok(graph) => graph,
        Err(error) if context.skip_unparseable => {
            return Ok(RowComparison::SkippedParse {
                row: global_row,
                smiles: smiles_text.clone(),
                error: error.to_string(),
            });
        }
        Err(error) => {
            return Err(
                format!("failed to prepare row {global_row}: {smiles_text}: {error}").into(),
            );
        }
    };
    let mut mismatches = Vec::new();

    for case in context.bit_cases {
        let expected = shard.bit_row(case.key(), row)?;
        let observed = case.compute(&smiles, &graph, context.maccs)?;
        if expected != observed.as_slice() {
            mismatches.push(RowMismatch {
                key: case.key().to_owned(),
                row: global_row,
                smiles: smiles_text.clone(),
                expected: format_entries(&active_bits_from_packed(expected, case.fp_size())),
                observed: format_entries(&active_bits_from_packed(&observed, case.fp_size())),
            });
        }
    }

    for case in context.count_cases {
        let expected = shard.count_row(case.key(), row)?;
        let observed = case.compute(&graph);
        if expected != observed {
            mismatches.push(RowMismatch {
                key: case.key().to_owned(),
                row: global_row,
                smiles: smiles_text.clone(),
                expected: format_entries(&expected),
                observed: format_entries(&observed),
            });
        }
    }

    Ok(RowComparison::Compared(mismatches))
}

fn print_summary(summaries: &BTreeMap<String, CaseSummary>) {
    for (key, summary) in summaries {
        println!(
            "{key}: checked={} mismatches={}",
            summary.checked, summary.mismatches
        );
        for sample in &summary.samples {
            println!(
                "  row={} smiles={} expected={} observed={}",
                sample.row, sample.smiles, sample.expected, sample.observed
            );
        }
    }
}

fn print_skipped_summary(skipped_rows: usize, skipped_samples: &[SkippedSample]) {
    println!("skipped_unparseable_rows={skipped_rows}");
    for sample in skipped_samples {
        println!(
            "  row={} smiles={} error={}",
            sample.row, sample.smiles, sample.error
        );
    }
}

fn format_entries(entries: &[(usize, u32)]) -> String {
    let parts = entries
        .iter()
        .map(|(index, count)| {
            if *count == 1 {
                index.to_string()
            } else {
                format!("{index}:{count}")
            }
        })
        .collect::<Vec<_>>();
    format!("[{}]", parts.join(","))
}

fn read_npy_u8<R: Read + std::io::Seek>(
    archive: &mut ZipArchive<R>,
    name: &str,
) -> Result<U8Array, DynError> {
    let npy = read_npy(archive, name)?;
    if npy.descr != "|u1" {
        return Err(format!("{name} has dtype {}, expected |u1", npy.descr).into());
    }
    Ok(U8Array {
        shape: npy.shape,
        data: npy.data,
    })
}

fn read_npy_u32<R: Read + std::io::Seek>(
    archive: &mut ZipArchive<R>,
    name: &str,
) -> Result<Vec<u32>, DynError> {
    let npy = read_npy(archive, name)?;
    if npy.descr != "<u4" {
        return Err(format!("{name} has dtype {}, expected <u4", npy.descr).into());
    }
    Ok(npy
        .data
        .chunks_exact(4)
        .map(|chunk| u32::from_le_bytes(chunk.try_into().expect("chunk size is 4")))
        .collect())
}

fn read_npy_u64<R: Read + std::io::Seek>(
    archive: &mut ZipArchive<R>,
    name: &str,
) -> Result<Vec<u64>, DynError> {
    let npy = read_npy(archive, name)?;
    if npy.descr != "<u8" {
        return Err(format!("{name} has dtype {}, expected <u8", npy.descr).into());
    }
    Ok(npy
        .data
        .chunks_exact(8)
        .map(|chunk| u64::from_le_bytes(chunk.try_into().expect("chunk size is 8")))
        .collect())
}

#[derive(Debug)]
struct NpyArray {
    descr: String,
    shape: Vec<usize>,
    data: Vec<u8>,
}

fn read_npy<R: Read + std::io::Seek>(
    archive: &mut ZipArchive<R>,
    name: &str,
) -> Result<NpyArray, DynError> {
    let mut file = archive.by_name(name)?;
    let mut bytes = Vec::with_capacity(file.size() as usize);
    file.read_to_end(&mut bytes)?;

    let mut cursor = Cursor::new(bytes);
    let mut magic = [0_u8; 6];
    cursor.read_exact(&mut magic)?;
    if &magic != b"\x93NUMPY" {
        return Err(format!("{name} is not an NPY file").into());
    }

    let mut version = [0_u8; 2];
    cursor.read_exact(&mut version)?;
    let header_len = match version {
        [1, 0] => {
            let mut len = [0_u8; 2];
            cursor.read_exact(&mut len)?;
            u16::from_le_bytes(len) as usize
        }
        [2 | 3, 0] => {
            let mut len = [0_u8; 4];
            cursor.read_exact(&mut len)?;
            u32::from_le_bytes(len) as usize
        }
        _ => return Err(format!("{name} uses unsupported NPY version {version:?}").into()),
    };

    let mut header = vec![0_u8; header_len];
    cursor.read_exact(&mut header)?;
    let header = String::from_utf8(header)?;
    let descr = parse_header_string(&header, "'descr': ")?;
    let fortran_order = parse_header_bool(&header, "'fortran_order': ")?;
    if fortran_order {
        return Err(format!("{name} uses Fortran order arrays").into());
    }
    let shape = parse_header_shape(&header)?;
    let item_size = dtype_item_size(&descr)?;
    let expected_len = shape
        .iter()
        .try_fold(item_size, |acc, dimension| acc.checked_mul(*dimension))
        .ok_or_else(|| format!("{name} shape overflows usize"))?;

    let mut data = Vec::new();
    cursor.read_to_end(&mut data)?;
    if data.len() != expected_len {
        return Err(format!(
            "{name} data length mismatch: got {}, expected {expected_len}",
            data.len()
        )
        .into());
    }

    Ok(NpyArray { descr, shape, data })
}

fn parse_header_string(header: &str, key: &str) -> Result<String, DynError> {
    let start = header
        .find(key)
        .ok_or_else(|| format!("missing NPY header key {key}"))?
        + key.len();
    let quote = header.as_bytes()[start] as char;
    if quote != '\'' && quote != '"' {
        return Err(format!("NPY header key {key} is not followed by a quoted string").into());
    }
    let value_start = start + 1;
    let value_end = header[value_start..]
        .find(quote)
        .ok_or_else(|| format!("unterminated NPY header string for {key}"))?
        + value_start;
    Ok(header[value_start..value_end].to_owned())
}

fn parse_header_bool(header: &str, key: &str) -> Result<bool, DynError> {
    let start = header
        .find(key)
        .ok_or_else(|| format!("missing NPY header key {key}"))?
        + key.len();
    if header[start..].starts_with("False") {
        Ok(false)
    } else if header[start..].starts_with("True") {
        Ok(true)
    } else {
        Err(format!("NPY header key {key} is not a bool").into())
    }
}

fn parse_header_shape(header: &str) -> Result<Vec<usize>, DynError> {
    let key = "'shape': (";
    let start = header.find(key).ok_or("missing NPY header shape")? + key.len();
    let end = header[start..]
        .find(')')
        .ok_or("unterminated NPY header shape")?
        + start;
    let body = header[start..end].trim();
    if body.is_empty() {
        return Ok(Vec::new());
    }
    body.split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(str::parse)
        .collect::<Result<Vec<_>, _>>()
        .map_err(Into::into)
}

fn dtype_item_size(descr: &str) -> Result<usize, DynError> {
    match descr {
        "|u1" => Ok(1),
        "<u4" => Ok(4),
        "<u8" => Ok(8),
        other => Err(format!("unsupported NPY dtype {other}").into()),
    }
}

fn decode_smiles(offsets: &[u64], utf8: &[u8]) -> Result<Vec<String>, DynError> {
    let mut smiles = Vec::with_capacity(offsets.len().saturating_sub(1));
    for window in offsets.windows(2) {
        let start = usize::try_from(window[0])?;
        let end = usize::try_from(window[1])?;
        smiles.push(std::str::from_utf8(&utf8[start..end])?.to_owned());
    }
    Ok(smiles)
}

fn set_packed_bits(packed: &mut [u8], bits: impl Iterator<Item = usize>) {
    for bit in bits {
        packed[bit / 8] |= 1 << (bit % 8);
    }
}

fn active_bits_from_packed(packed: &[u8], bit_count: usize) -> Vec<(usize, u32)> {
    (0..bit_count)
        .filter(|bit| packed[*bit / 8] & (1 << (*bit % 8)) != 0)
        .map(|bit| (bit, 1))
        .collect()
}

fn byte_width(bit_count: usize) -> usize {
    bit_count.div_ceil(8)
}

impl fmt::Display for MismatchSample {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "row={} smiles={} expected={} observed={}",
            self.row, self.smiles, self.expected, self.observed
        )
    }
}
