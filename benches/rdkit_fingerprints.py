from __future__ import annotations

import argparse
import gzip
import statistics
import time
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdFingerprintGenerator as FG

RDLogger.DisableLog("rdApp.*")

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CORPUS = REPO_ROOT / "tests" / "fixtures" / "pubchem_benchmark_10000_smiles.txt.gz"

ECFP_CASES = [(0, 64), (1, 128), (2, 128), (2, 2048), (4, 2048), (5, 2048), (5, 4096)]
FOLDED_BIT_SIZES = [64, 128, 256, 512, 1024, 2048, 4096]


def load_smiles(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def load_molecules(path: Path) -> list[Any]:
    smiles_list = load_smiles(path)
    molecules = []
    for smiles in smiles_list:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise RuntimeError(f"RDKit failed to parse benchmark SMILES: {smiles}")
        molecules.append(molecule)
    return molecules


def morgan_bit_case(radius: int, fp_size: int) -> tuple[str, Callable[[Any], Any], Callable[[Any], int]]:
    generator = FG.GetMorganGenerator(
        radius=radius,
        fpSize=fp_size,
        includeChirality=False,
        useBondTypes=True,
        includeRingMembership=True,
    )
    return (
        f"ecfp_bits/r{radius}_n{fp_size}",
        generator.GetFingerprint,
        lambda fingerprint: fingerprint.GetNumOnBits(),
    )


def morgan_count_case(radius: int, fp_size: int) -> tuple[str, Callable[[Any], Any], Callable[[Any], int]]:
    generator = FG.GetMorganGenerator(
        radius=radius,
        fpSize=fp_size,
        includeChirality=False,
        useBondTypes=True,
        includeRingMembership=True,
    )
    return (
        f"ecfp_counts/r{radius}_n{fp_size}",
        generator.GetCountFingerprint,
        lambda fingerprint: sum(fingerprint.GetNonzeroElements().values()),
    )


def atom_pair_case(fp_size: int) -> tuple[str, Callable[[Any], Any], Callable[[Any], int]]:
    generator = FG.GetAtomPairGenerator(
        minDistance=1,
        maxDistance=30,
        includeChirality=False,
        use2D=True,
        countSimulation=True,
        fpSize=fp_size,
    )
    return (
        f"atom_pair/n{fp_size}",
        generator.GetFingerprint,
        lambda fingerprint: fingerprint.GetNumOnBits(),
    )


def topological_torsion_case(fp_size: int) -> tuple[str, Callable[[Any], Any], Callable[[Any], int]]:
    generator = FG.GetTopologicalTorsionGenerator(
        includeChirality=False,
        torsionAtomCount=4,
        countSimulation=True,
        fpSize=fp_size,
    )
    return (
        f"topological_torsion/n{fp_size}",
        generator.GetFingerprint,
        lambda fingerprint: fingerprint.GetNumOnBits(),
    )


def rdk_case(fp_size: int) -> tuple[str, Callable[[Any], Any], Callable[[Any], int]]:
    def compute(molecule: Any) -> Any:
        return Chem.RDKFingerprint(
            molecule,
            minPath=1,
            maxPath=7,
            fpSize=fp_size,
            nBitsPerHash=2,
            useHs=True,
            tgtDensity=0.0,
            minSize=128,
            branchedPaths=True,
            useBondOrder=True,
        )

    return (f"rdk/n{fp_size}", compute, lambda fingerprint: fingerprint.GetNumOnBits())


def maccs_case() -> tuple[str, Callable[[Any], Any], Callable[[Any], int]]:
    return ("maccs/n167", MACCSkeys.GenMACCSKeys, lambda fingerprint: fingerprint.GetNumOnBits())


def selected_cases(
    fingerprint: str,
) -> Iterable[tuple[str, Callable[[Any], Any], Callable[[Any], int]]]:
    if fingerprint in {"all", "ecfp"}:
        for radius, fp_size in ECFP_CASES:
            yield morgan_bit_case(radius, fp_size)
    if fingerprint in {"all", "counted-ecfp"}:
        for radius, fp_size in ECFP_CASES:
            yield morgan_count_case(radius, fp_size)
    if fingerprint in {"all", "atom-pair"}:
        for fp_size in FOLDED_BIT_SIZES:
            yield atom_pair_case(fp_size)
    if fingerprint in {"all", "topological-torsion"}:
        for fp_size in FOLDED_BIT_SIZES:
            yield topological_torsion_case(fp_size)
    if fingerprint in {"all", "rdk"}:
        for fp_size in FOLDED_BIT_SIZES:
            yield rdk_case(fp_size)
    if fingerprint in {"all", "maccs"}:
        yield maccs_case()


def benchmark_case(
    molecules: list[Any],
    compute: Callable[[Any], Any],
    checksum: Callable[[Any], int],
    warmup: int,
    repeats: int,
) -> tuple[float, float, int]:
    for _ in range(warmup):
        for molecule in molecules:
            checksum(compute(molecule))

    durations_ns = []
    checksum_value = 0
    for _ in range(repeats):
        started = time.perf_counter_ns()
        for molecule in molecules:
            checksum_value += checksum(compute(molecule))
        durations_ns.append(time.perf_counter_ns() - started)

    median_ns = statistics.median(durations_ns)
    mean_ns = statistics.fmean(durations_ns)
    return median_ns / 1_000_000.0, mean_ns / 1_000_000.0, checksum_value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fingerprint",
        choices=("all", "ecfp", "counted-ecfp", "atom-pair", "topological-torsion", "rdk", "maccs"),
        default="all",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--warmup", type=int, default=3)
    parser.add_argument("--repeats", type=int, default=10)
    args = parser.parse_args()

    molecules = load_molecules(args.input)
    print(f"RDKit benchmark corpus: {len(molecules)} molecules")
    print(f"input={args.input}")
    print(f"warmup={args.warmup} repeats={args.repeats}")
    print("case,median_ms_per_corpus,mean_ms_per_corpus,median_us_per_molecule,checksum")

    for name, compute, checksum in selected_cases(args.fingerprint):
        median_ms, mean_ms, checksum_value = benchmark_case(
            molecules,
            compute,
            checksum,
            args.warmup,
            args.repeats,
        )
        print(
            f"{name},{median_ms:.3f},{mean_ms:.3f},"
            f"{median_ms / len(molecules) * 1000.0:.3f},{checksum_value}"
        )


if __name__ == "__main__":
    main()
