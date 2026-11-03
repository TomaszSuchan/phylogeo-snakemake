#!/usr/bin/env bash
# Run OrientAGraph, retrying with a bumped seed on the known outgroup re-root crash:
#   ERROR: Calling get_child_node_mig on a non-migration node
#
# Usage:
#   run_orientagraph_retry.sh LOG PREFIX FRQ M K SEED MAX_TRIES [extra orientagraph flags...]
set -euo pipefail

log=$1
prefix=$2
frq=$3
m=$4
k=$5
seed=$6
max_tries=$7
shift 7

mkdir -p "$(dirname "$log")" "$(dirname "$prefix")"
: > "$log"

attempt=1
while (( attempt <= max_tries )); do
  rm -f \
    "${prefix}.llik" \
    "${prefix}.treeout.gz" \
    "${prefix}.cov.gz" \
    "${prefix}.covse.gz" \
    "${prefix}.modelcov.gz" \
    "${prefix}.vertices.gz" \
    "${prefix}.edges.gz" \
    "${prefix}.mltree"

  {
    echo "=== OrientAGraph attempt ${attempt}/${max_tries} seed=${seed} ==="
    echo "command: orientagraph -i ${frq} -o ${prefix} -m ${m} -k ${k} -seed ${seed} $*"
  } >> "$log"

  set +e
  orientagraph -i "$frq" -o "$prefix" -m "$m" -k "$k" -seed "$seed" "$@" >> "$log" 2>&1
  status=$?
  set -e

  if (( status == 0 )); then
    exit 0
  fi

  if grep -q 'Calling get_child_node_mig on a non-migration node' "$log"; then
    echo "OrientAGraph outgroup re-root crash at seed=${seed}; retrying with a new seed" >> "$log"
    seed=$((seed + 7919 * attempt))
    attempt=$((attempt + 1))
    continue
  fi

  echo "OrientAGraph failed with exit status ${status} (not a known re-root crash)" >> "$log"
  exit "$status"
done

echo "OrientAGraph gave up after ${max_tries} seeds (repeated get_child_node_mig crash)" >> "$log"
exit 1
