#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/ams_run_common.sh"

# Matches the 20251228 R-MAT panels: ER 19-21, G500 17-19; edge factors 2-32.
# The plotted run uses NUMA mask 1. Its AMS setup used 64 logical CPUs local
# to node 0, with ten iterations including the warmup.
RUN_NUM=${RUN_NUM:-10}
THREAD_NUM=${THREAD_NUM:-64}
AMS_HIGH_THREAD_OVERRIDE=${AMS_HIGH_THREAD_NUM:-}
CORE_SET=${CORE_SET:-0-31,64-95}
LOG_DIR=${LOG_DIR:-"$ROOT/logs/ams-rmat"}
mkdir -p "$LOG_DIR"

for family in er g500; do
	for bin in "${AMS_BINS[@]}"; do
		: > "$LOG_DIR/${family}_1_${THREAD_NUM}t_${bin}.log"
	done
done


for bin in "${AMS_BINS[@]}"; do
	case "$bin" in
		ams_hashspgemm)
			high_threads=16; prt=0.90; acc=5.00; smt_symb=11000000; smt_numc=2500000
			;;
		ams_hserscspgemm)
			high_threads=22; prt=0.85; acc=5.00; smt_symb=7000000; smt_numc=2100000
			;;
		ams_hybaccspgemm)
			high_threads=22; prt=0.85; acc=2.00; smt_symb=7000000; smt_numc=2100000
			;;
	esac
	if [[ -n "$AMS_HIGH_THREAD_OVERRIDE" ]]; then high_threads=$AMS_HIGH_THREAD_OVERRIDE; fi
	build_ams_targets \
		"$THREAD_NUM" ON 1 OFF config/mc/mxc0.cfg "$high_threads" \
		"$prt" "$acc" "$smt_symb" "$smt_numc" "$RUN_NUM" "$bin"

	for family in er g500; do
		if [[ "$family" == er ]]; then scales=(19 20 21); else scales=(17 18 19); fi
		for scale in "${scales[@]}"; do
			for edge in 2 4 8 16 32; do
				run_ams_dataset "${family}_${scale}_${edge}" "$CORE_SET" \
					"${family}_1_${THREAD_NUM}t" "$bin"
			done
		done
	done
done

echo "R-MAT experiment logs: $LOG_DIR"
