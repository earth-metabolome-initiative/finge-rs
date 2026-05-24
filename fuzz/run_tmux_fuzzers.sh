#!/usr/bin/env bash
set -euo pipefail

SESSION_NAME="finge-fuzz"
ATTACH="${FUZZ_ATTACH:-1}"

while (($#)); do
  case "$1" in
    --replace)
      # Kept for backwards compatibility — replacement is now the default.
      shift
      ;;
    --no-attach)
      ATTACH=0
      shift
      ;;
    *)
      SESSION_NAME="$1"
      shift
      ;;
  esac
done

if ! command -v tmux >/dev/null 2>&1; then
  echo "tmux is required but was not found in PATH" >&2
  exit 1
fi

if ! command -v cargo >/dev/null 2>&1; then
  echo "cargo is required but was not found in PATH" >&2
  exit 1
fi

if ! cargo fuzz --help >/dev/null 2>&1; then
  cat >&2 <<'EOF'
cargo-fuzz is required but does not appear to be installed.

Install it with:
  cargo install cargo-fuzz
EOF
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# libFuzzer runtime knobs (passed after `--` to cargo-fuzz). Mirror the
# subql `fuzz.sh` defaults; finge-rs-specific rationale:
#   -timeout=15        abort a single input after 15s. The smiles-parser
#                      `perceive_aromaticity_for(RdkitDefault)` path is
#                      pathologically slow (~16s) on certain adversarial
#                      69-atom SMILES with reused small-digit ring closures
#                      across many disconnected fragments. Without this,
#                      a single slow input lets libFuzzer's RSS drift past
#                      the default 2 GiB rss_limit and trip code 71.
#                      (Tracked upstream in smiles-parser; until fixed,
#                      the timeout converts these from OOMs into clean
#                      timeout artifacts that the fuzzer can skip past.)
#   -max_len=512       cap generated input size. `parse_smiles` already
#                      rejects strings longer than 128 bytes, so a 512-byte
#                      libfuzzer cap leaves ample room for the `arbitrary`
#                      String decoding overhead without wasting cycles on
#                      multi-KiB byte streams that the harness drops.
#   -rss_limit_mb=8192 raise libFuzzer's RSS ceiling from 2 GiB to 8 GiB.
#                      libFuzzer accumulates corpus / coverage state across
#                      iterations; on this machine 8 GiB is ample headroom
#                      and prevents code-71 false alarms when an adversarial
#                      input briefly spikes RSS during a long session.
LIBFUZZER_DEFAULTS="-timeout=15 -max_len=512 -rss_limit_mb=8192"
COMMON_ARGS="${LIBFUZZER_DEFAULTS} ${FUZZ_RUN_ARGS:-}"

if tmux has-session -t "${SESSION_NAME}" 2>/dev/null; then
  echo "killing existing tmux session '${SESSION_NAME}'" >&2
  tmux kill-session -t "${SESSION_NAME}"
fi

# Window 0: the four fingerprint-family fuzzers (one tile per family).
tmux new-session -d -s "${SESSION_NAME}" -n "fingerprints" -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run ecfp -- ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0" -h -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run atom_pair -- ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0.0" -v -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run topological_torsion -- ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0.1" -v -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run maccs -- ${COMMON_ARGS}'"

tmux select-layout -t "${SESSION_NAME}:0" tiled

# Window 1: the eight per-mutator invariance fuzzers plus the
# MutatorMix composition fuzzer (one tile per target).
tmux new-window -t "${SESSION_NAME}" -n "mutators" -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run mutator_atomic_number -- ${COMMON_ARGS}'"

for target in mutator_hypervalent mutator_h_count mutator_formal_charge \
              mutator_isotope mutator_ring_flag mutator_bond_type \
              mutator_topology mutator_mix; do
  tmux split-window -t "${SESSION_NAME}:1" -c "${REPO_ROOT}" \
    "bash -lc 'cargo fuzz run ${target} -- ${COMMON_ARGS}'"
  tmux select-layout -t "${SESSION_NAME}:1" tiled
done

tmux select-window -t "${SESSION_NAME}:0"
tmux set-option -t "${SESSION_NAME}" remain-on-exit on

cat <<EOF
Started tmux session '${SESSION_NAME}' with two windows:
  window 0 'fingerprints' (4 panes):
    - ecfp
    - atom_pair
    - topological_torsion
    - maccs
  window 1 'mutators' (9 panes):
    - mutator_atomic_number
    - mutator_hypervalent
    - mutator_h_count
    - mutator_formal_charge
    - mutator_isotope
    - mutator_ring_flag
    - mutator_bond_type
    - mutator_topology
    - mutator_mix

Re-attach later with:
  tmux attach -t ${SESSION_NAME}

Switch windows: prefix + n (next) / prefix + p (prev) / prefix + 0/1

libFuzzer defaults (appended after \`--\` to every \`cargo fuzz run\`):
  ${LIBFUZZER_DEFAULTS}

Optional:
  FUZZ_RUN_ARGS='-max_total_time=300' ${BASH_SOURCE[0]} ${SESSION_NAME}
    # appended after the defaults, so it can both override and extend them
  FUZZ_ATTACH=0 ${BASH_SOURCE[0]} ${SESSION_NAME}    # start without attaching
  ${BASH_SOURCE[0]} --no-attach ${SESSION_NAME}      # same, via flag
EOF

if [[ "${ATTACH}" == "1" ]]; then
  if [[ -t 1 ]]; then
    exec tmux attach -t "${SESSION_NAME}"
  else
    echo "stdout is not a TTY; skipping auto-attach" >&2
    echo "Attach manually with: tmux attach -t ${SESSION_NAME}" >&2
  fi
fi
