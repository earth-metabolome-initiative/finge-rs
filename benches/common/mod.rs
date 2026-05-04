use std::{
    env,
    fs::File,
    io::Read,
    path::{Path, PathBuf},
    time::Duration,
};

use criterion::{BenchmarkGroup, Throughput, measurement::WallTime};
use flate2::read::GzDecoder;
use smiles_parser::smiles::Smiles;

const DEFAULT_BENCH_CORPUS: &str = "tests/fixtures/pubchem_benchmark_10000_smiles.txt.gz";

pub fn load_smiles_corpus() -> Vec<Smiles> {
    let path = env::var_os("FINGE_RS_BENCH_SMILES")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(DEFAULT_BENCH_CORPUS));
    let mut contents = String::new();
    open_corpus(&path)
        .read_to_string(&mut contents)
        .unwrap_or_else(|error| {
            panic!(
                "failed to read benchmark corpus {}: {error}",
                path.display()
            )
        });

    contents
        .lines()
        .filter(|line| !line.is_empty())
        .map(|smiles| {
            smiles.parse::<Smiles>().unwrap_or_else(|error| {
                panic!("benchmark corpus SMILES should parse: {smiles}: {error}")
            })
        })
        .collect()
}

fn open_corpus(path: &Path) -> Box<dyn Read> {
    let file = File::open(path).unwrap_or_else(|error| {
        panic!(
            "failed to open benchmark corpus {}: {error}",
            path.display()
        )
    });
    if path
        .extension()
        .is_some_and(|extension| extension == std::ffi::OsStr::new("gz"))
    {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    }
}

pub fn configure_group(group: &mut BenchmarkGroup<'_, WallTime>, corpus_len: usize) {
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(4));
    group.warm_up_time(Duration::from_secs(1));
    group.throughput(Throughput::Elements(corpus_len as u64));
}
