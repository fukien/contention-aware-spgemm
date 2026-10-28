#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/ams_run_common.sh"

# Run all three datasets with the mask 1 profile: local-node placement and no
# weighted interleaving. The runner uses 64 logical CPUs and ten iterations.
RUN_NUM=${RUN_NUM:-10}
THREAD_NUM=${THREAD_NUM:-64}
AMS_HIGH_THREAD_NUM=${AMS_HIGH_THREAD_NUM:-22}
LOG_DIR=${LOG_DIR:-"$ROOT/logs/ams-join-aggregate"}
mkdir -p "$LOG_DIR"

for dataset in imdb_movie bkgn chbc_s_uni; do
	for bin in "${AMS_BINS[@]}"; do
		: > "$LOG_DIR/${dataset}_${bin}.log"
	done
done

# All datasets: local-node placement, no interleaving, mask 1.
build_ams_targets \
	"$THREAD_NUM" ON 1 OFF config/mc/mxc0.cfg "$AMS_HIGH_THREAD_NUM" \
	0.85 1.2116 7000000 2100000 "$RUN_NUM"
run_ams_dataset imdb_movie 0-31,64-95 imdb_movie
run_ams_dataset bkgn 0-31,64-95 bkgn
run_ams_dataset chbc_s_uni 0-31,64-95 chbc_s_uni

echo "Join-Aggregate experiment logs: $LOG_DIR"
