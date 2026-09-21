#!/bin/bash
# Verify cellranger_localcores/localmem derive from the Slurm allocation rather
# than from hardcoded numbers. The pre-restructure scripts paired
# "#SBATCH --mem=12G" with "--localmem=160" and "--mem=10G" with
# "--localmem 100", both of which get the job OOM-killed.
#
# Usage: bash tests/test_config_memory.sh
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
RESULTS=$(mktemp)
trap 'rm -f "$RESULTS"' EXIT
# Counters are written to a file because each scenario runs in a subshell to
# isolate its exported SLURM_* variables.
export RESULTS
check() { # desc expected actual
  if [[ "$2" == "$3" ]]; then echo "  PASS  $1 -> $3"; echo pass >> "$RESULTS"
  else echo "  FAIL  $1: expected $2, got $3"; echo fail >> "$RESULTS"; fi
}
note() { echo "  $1  $2"; echo "$3" >> "$RESULTS"; }
export -f check note

echo "=== 125G / 16 cpus (02_count.sh allocation) ==="
( export SLURM_CPUS_PER_TASK=16 SLURM_MEM_PER_NODE=128000
  source ./config.sh
  check "localcores" 16 "$(cellranger_localcores)"
  check "localmem (128000MB -> 125G - 4 headroom)" 121 "$(cellranger_localmem)" )

echo "=== 64G / 16 cpus (01_mkfastq.sh allocation) ==="
( export SLURM_CPUS_PER_TASK=16 SLURM_MEM_PER_NODE=65536
  source ./config.sh
  check "localmem" 60 "$(cellranger_localmem)" )

echo "=== mem-per-cpu instead of mem-per-node ==="
( export SLURM_CPUS_PER_TASK=8 SLURM_MEM_PER_CPU=4096
  unset SLURM_MEM_PER_NODE
  source ./config.sh
  check "localcores" 8 "$(cellranger_localcores)"
  check "localmem (8 x 4096MB = 32G - 4)" 28 "$(cellranger_localmem)" )

echo "=== the original bug: localmem can no longer exceed the allocation ==="
( export SLURM_CPUS_PER_TASK=16 SLURM_MEM_PER_NODE=12288   # the old --mem=12G
  source ./config.sh
  got=$(cellranger_localmem)
  # Old script passed --localmem=160 against this 12G allocation.
  if (( got < 12 )); then note PASS "12G allocation yields ${got}G, not 160G" pass
  else note FAIL "got ${got}G from a 12G allocation" fail; fi )

echo "=== outside Slurm: refuses to guess ==="
( unset SLURM_MEM_PER_NODE SLURM_MEM_PER_CPU SLURM_CPUS_PER_TASK SLURM_CPUS_ON_NODE
  source ./config.sh
  if cellranger_localmem >/dev/null 2>&1; then note FAIL "guessed a memory value" fail
  else note PASS "errored rather than guessing" pass; fi )

echo "=== allocation too small for headroom ==="
( export SLURM_MEM_PER_NODE=2048
  source ./config.sh
  if cellranger_localmem >/dev/null 2>&1; then note FAIL "returned a value" fail
  else note PASS "errored on a 2G allocation" pass; fi )

pass=$(grep -c '^pass$' "$RESULTS" || true)
fail=$(grep -c '^fail$' "$RESULTS" || true)
echo ""
echo "passed: $pass   failed: $fail"
[[ "$fail" -eq 0 ]]
