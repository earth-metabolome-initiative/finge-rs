# smiles-parser: `perceive_aromaticity_for(RdkitDefault)` pathological slowdown

**Status: FIXED upstream on 2026-05-24 at revision `ba32bd63`.** The original
report and post-fix timings are kept below for reference.

Discovered via finge-rs fuzzing on 2026-05-22, profiled and written up
2026-05-24. Original timings tested against `smiles-parser` revision
`0dbb557e` (the bracket-atom valence fix from 2026-05-23 with the
aromaticity perf bug still present).

## Bug

`Smiles::perceive_aromaticity_for(AromaticityPolicy::RdkitDefault)` exhibits
superlinear (likely exponential in cycle count) runtime on small adversarial
inputs. Release-build timings on fuzz-discovered samples sit between
**840 ms and 5.6 seconds for 58 to 89 atoms**. Debug-build timings on the
same inputs reach **16 to 26 seconds**. Downstream operations
(`kekulize_standalone`, `to_invariant`, fingerprinting) run in microseconds
on the same graphs.

## Severity

Any fuzz harness or untrusted-input pipeline that goes through
`Smiles::try_prepare` (which internally calls `perceive_aromaticity_for(RdkitDefault)`)
becomes practically unfuzzable when adversarial inputs surface. Without
per-input timeouts the fuzz session stalls. With timeouts it accumulates
`timeout-*` artifacts at a rate that drowns out real signal.

## Minimal reproducer

Drop into a release-mode binary or a `#[test]` run with `--release`:

```rust
use smiles_parser::smiles::{AromaticityPolicy, Smiles};
use smiles_parser::prelude::MolecularGraph;

let smiles: Smiles = "bFFFFFFFFn3SS752sBS5S73B1S3=5S7B2S5S73B2(cS3n5S7B2S5S732BS3n5S7B2S5S73B27S53SnB1S5S7B33Sn2)5S7B2C5S73S53BB27n2S5S7F3-2"
    .parse()
    .unwrap();
let base = smiles.kekulize_standalone().unwrap_or_else(|_| smiles.clone());
let _ = base.perceive_aromaticity_for(AromaticityPolicy::RdkitDefault);
// ~5.6 seconds release-build on 58 atoms.
```

## Profile table

Release build (`cargo test --lib --release`), all on a Threadripper
5975WX. `len` is SMILES bytes, `atoms` is post-parse atom count, `frags`
is the number of `.`-separated fragments.

| len | atoms | frags | kekulize | perceive_aromaticity | SMILES |
|---:|---:|---:|---:|---:|---|
| 118 | 58 | 1 | 11 us | **5.577 s** | `bFFFFFFFFn3SS752sBS5S73B1S3=5S7B2S5S73B2(cS3n5S7B2S5S732BS3n5S7B2S5S73B27S53SnB1S5S7B33Sn2)5S7B2C5S73S53BB27n2S5S7F3-2` |
| 127 | 85 | 1 | 13 us | **3.040 s** | `C2CFCBBC1FB2B2nBBBFCF22BBBnF$CF11BBBI-FFB2B2nCB2B2nBBB22BBFB2B211BBC1F2B2nCBF2B2FF2B2nF2B2nCBF2B2FF2B2nFB2B2nB2B2FB2B2nB2B2FBC2` |
| 104 | 61 | 1 | 11 us | **1.778 s** | `F3oFF12s3FP3PS2SFF11P2F1s2SP2F1s2S1F/sSFF12s3FP3PS2SFF11P2F1s2SP2F1s2S1F/s1FS2FS1FP2Fs1FS2FS1F1SbS3#SS21` |
| 110 | 73 | 1 | 12 us | **1.343 s** | `F789C/C-C3C7$Cn3CC7CC3nCC7nC3CF7CC3nCC7nC30C7CC3nCC7nC3CF7CC3nCCNNF8CNBCS7nC3C7CC3nCC7nC3CF7CC3nC0C7CC39CCnNCN` |
| 126 | 69 | 2 | 28 us | **1.282 s** | `C7pC02CB1CFB2N211Bs.CF2N211CF22BF/FB12N12BBsC11BBCF22BF/C0BB7c2B1B7PcF#pCBBC210CF2N11FC2BB2N211BBCFCFB2N211/C0BB7c2B1B2B[Os]=2` |
| 121 | 72 | 4 | 12 us | **850 ms** | `FF9p7FFcFFFnF8#F7pF88FF.sF#F7pF8FsBB8#F7pF8Fs8FF#F7pFBB81[Ac]F.F8#F7pF88FF.sF#F7pF8FsBB8#N7pF8Fs8FF#F7pB=81F77Cp8FF9FcF87` |
| 128 | 89 | 1 | 20 us | **843 ms** | `S8IoNNIII4II5oI88I4IIIINNI5II5NI88IIIIIIII5NNI5NI88IINNII84NNN5=I5I8IoNI8No8N4I4IIINN5I5IINNN88IoI4II5oI88O4IINNI5IINNI48N5I8NI8` |

Atom count alone is a poor predictor — the 58-atom input is 6x slower
than the 89-atom one. The slow inputs share structural features that
appear to be the real driver.

## Surface characteristics of the trigger

Every slow input has at least three of the following:

- Heavy reuse of small-digit ring closures (`1`, `2`, `3`, `7`, `8`, `0`)
  across positions, often with the same digit closing multiple distinct
  rings within one fragment.
- Mixed aromatic-letter atoms (`p`, `c`, `b`, `s`, `n`, `o`) that would
  not normally be aromatic outside specific cycles, scattered through
  aliphatic context.
- Adjacent `B` runs (boron chains) creating dense low-electronegativity
  cycle candidates.
- Bond-type punctuation (`=`, `#`, `/`, `\`, `$`, `-`) sprinkled
  mid-fragment, opening many alternative kekule assignments.
- A bracket atom with an unusual element and a ring closure attached
  (`[Os]=2`, `[Ac]F`, `[49S]0`, `[205Pb]C#`, `[Sc]0`) sitting near other
  ring-closure digits.

The combination produces graphs with many small candidate cycles where
each cycle's aromaticity assignment depends on its neighbours'
assignments. The timings (1.3 s for 2 fragments at 69 atoms versus 5.6 s
for 1 fragment at 58 atoms) suggest the cost is dominated by per-fragment
cycle enumeration or per-cycle backtracking rather than by total atom
count.

## Suggested fix direction

Profile a single representative input (the 58-atom 5.6-second case is the
cleanest) with `perf` or `flamegraph` against a release build to identify
which subroutine inside `perceive_aromaticity_for(RdkitDefault)` dominates.
Candidates to look at first:

1. SSSR (smallest set of smallest rings) or minimum-cycle-basis
   enumeration. For graphs with many shared edges across cycles, naive
   enumeration can blow up combinatorially.
2. The aromatic-ring assignment fixpoint. If the algorithm iterates "for
   each ring, decide aromaticity given neighbours' current state" to a
   fixpoint, pathological topologies can keep flipping assignments.
3. Per-ring backtracking through kekule alternatives. The `=`, `#`, `$`
   punctuation in the slow inputs interacts here.
4. Caching: if intermediate results (ring membership, pi-electron counts,
   etc.) are not memoised across the algorithm's passes, repeated
   recomputation may dominate.

A short-circuit fix that is likely acceptable for fuzz-driven use cases:
bail out with `Err(_)` when SSSR cardinality (or candidate-aromatic-ring
count) exceeds some threshold per fragment. finge-rs's `try_prepare`
already propagates the error and the fuzz harness handles it via `return`,
so `Err` is a safe outcome from the user's perspective. The threshold can
be generous (50 rings per fragment, say) and still catch every slow input
in this corpus while letting all real molecules through.

## Repro files (verbatim corpus)

In `finge-rs/fuzz/artifacts/`:

- `topological_torsion/oom-2785e9133fafd0cd34dd0a16b3f70df4810740e1`
- `topological_torsion/oom-d89e1c8d711221e18305adc7b172156df94c2249`
- `topological_torsion/oom-f6c5047f28c1f88f02be3e1f88f06e66ff93822f`
- `ecfp/timeout-50af30a1488d76e7ccdfde2ddf88f46233efca81`
- `mutator_h_count/timeout-32a342f8a61e07d8d9d514ef7e1777e3fb49fba0`
- `mutator_formal_charge/timeout-4536b72c30e4d81ae5f46dd4fdf6bbfad203c7eb`
- `mutator_isotope/timeout-90d95355f6a3adfa6147793b41cbee1ffd3687d1`

Each decodes to the corresponding SMILES via `cargo +nightly fuzz fmt <target> <artifact>`.
They can be lifted directly into a `smiles-parser` regression test.

## finge-rs side workaround already deployed

`-timeout=15` + `-timeout_exitcode=0` in `fuzz/run_tmux_fuzzers.sh`
converts these into clean `timeout-*` artifacts without stopping the fuzz
session. Will be removed once the upstream perf fix lands.

## Post-fix verification (2026-05-24, revision `ba32bd63`)

Re-running the same release-mode profile after the upstream fix:

| atoms | before | after | speedup |
|---:|---:|---:|---:|
| 58 | 5.577 s | 47 ms | **118x** |
| 85 | 3.040 s | 110 ms | 28x |
| 72 | 850 ms | 138 ms | 6x |
| 61 | 1.778 s | 91 ms | 20x |
| 89 | 843 ms | 12 ms | 70x |
| 73 | 1.343 s | 89 ms | 15x |
| 69 | 1.282 s | 25 ms | 51x |

Total corpus runtime: ~14 s → ~510 ms. Slowest case is 138 ms, well under
any reasonable per-input fuzz timeout.

Re-replaying every `timeout-*` and `oom-*` artifact in
`fuzz/artifacts/` through `cargo +nightly fuzz run` with the
launcher's `-timeout=15 -rss_limit_mb=8192` confirms all finish under 4 s
wall (mostly libFuzzer startup; per-input runtime is sub-second). None
re-trip the timeout or RSS thresholds.

The `-timeout=15` + `-timeout_exitcode=0` knobs in the launcher are kept
as defense-in-depth against future pathological inputs (same rationale
as subql's `fuzz.sh`), but the smiles-parser-specific reason that
motivated them is resolved.
