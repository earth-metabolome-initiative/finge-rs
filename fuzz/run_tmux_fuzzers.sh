#!/usr/bin/env bash
set -euo pipefail

SESSION_NAME="finge-fuzz"
REPLACE_SESSION="${FUZZ_REPLACE_SESSION:-0}"

while (($#)); do
  case "$1" in
    --replace)
      REPLACE_SESSION=1
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
  if [[ "${REPLACE_SESSION}" == "1" ]]; then
    tmux kill-session -t "${SESSION_NAME}"
  else
    pane_states="$(tmux list-panes -t "${SESSION_NAME}" -F '#{pane_dead}' 2>/dev/null || true)"
    if [[ -n "${pane_states}" ]] && ! grep -q '^0$' <<<"${pane_states}"; then
      echo "replacing stale tmux session '${SESSION_NAME}' with only dead panes" >&2
      tmux kill-session -t "${SESSION_NAME}"
    else
      echo "tmux session '${SESSION_NAME}' already exists" >&2
      echo "Attach with: tmux attach -t ${SESSION_NAME}" >&2
      echo "Or replace it with: ${BASH_SOURCE[0]} --replace ${SESSION_NAME}" >&2
      exit 1
    fi
  fi
fi

tmux new-session -d -s "${SESSION_NAME}" -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run ecfp ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0" -h -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run atom_pair ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0.0" -v -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run topological_torsion ${COMMON_ARGS}'"

tmux split-window -t "${SESSION_NAME}:0.1" -v -c "${REPO_ROOT}" \
  "bash -lc 'cargo fuzz run maccs ${COMMON_ARGS}'"

tmux select-layout -t "${SESSION_NAME}:0" tiled
tmux set-option -t "${SESSION_NAME}" remain-on-exit on

cat <<EOF
Started tmux session '${SESSION_NAME}' with 4 fuzzers:
  - ecfp
  - atom_pair
  - topological_torsion
  - maccs

Attach with:
  tmux attach -t ${SESSION_NAME}

Optional:
  FUZZ_RUN_ARGS='-max_total_time=300' ${BASH_SOURCE[0]} ${SESSION_NAME}
  ${BASH_SOURCE[0]} --replace ${SESSION_NAME}
EOF
