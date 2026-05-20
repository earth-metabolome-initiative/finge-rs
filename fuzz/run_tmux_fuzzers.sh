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
COMMON_ARGS="${FUZZ_RUN_ARGS:-}"

if tmux has-session -t "${SESSION_NAME}" 2>/dev/null; then
  echo "killing existing tmux session '${SESSION_NAME}'" >&2
  tmux kill-session -t "${SESSION_NAME}"
fi

# Window 0: the four fingerprint-family fuzzers (one tile per family).
tmux new-session -d -s "${SESSION_NAME}" -n "fingerprints" -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run ecfp ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0" -h -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run atom_pair ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0.0" -v -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run topological_torsion ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0.1" -v -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run maccs ${COMMON_ARGS}'"

tmux select-layout -t "${SESSION_NAME}:0" tiled

# Window 1: the eight per-mutator invariance fuzzers (one tile per mutator).
tmux new-window -t "${SESSION_NAME}" -n "mutators" -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run mutator_atomic_number ${COMMON_ARGS}'"

for target in mutator_hypervalent mutator_h_count mutator_formal_charge \
              mutator_isotope mutator_ring_flag mutator_bond_type mutator_topology; do
  tmux split-window -t "${SESSION_NAME}:1" -c "${REPO_ROOT}" \
    "bash -lc 'cargo fuzz run ${target} ${COMMON_ARGS}'"
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
  window 1 'mutators' (8 panes):
    - mutator_atomic_number
    - mutator_hypervalent
    - mutator_h_count
    - mutator_formal_charge
    - mutator_isotope
    - mutator_ring_flag
    - mutator_bond_type
    - mutator_topology

Re-attach later with:
  tmux attach -t ${SESSION_NAME}

Switch windows: prefix + n (next) / prefix + p (prev) / prefix + 0/1

Optional:
  FUZZ_RUN_ARGS='-max_total_time=300' ${BASH_SOURCE[0]} ${SESSION_NAME}
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
