#!/usr/bin/env bash

# Shared build/run helpers for the R-MAT and Join-Aggregate experiment scripts.
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
LOG_DIR=${LOG_DIR:-"$ROOT/logs"}
BUILD_JOBS=${BUILD_JOBS:-4}
NUMACTL=${NUMACTL:-numactl}
AMS_BINS=(ams_hashspgemm ams_hserscspgemm ams_hybaccspgemm)

build_ams_targets() {
	local threads=$1 hyper=$2 mask=$3 weighted=$4 cfg=$5 high_threads=$6
	local prt=$7 acc=$8 smt_symb=$9 smt_numc=${10} run_num=${11}
	local hyper_opt=OFF weighted_opt=OFF smt_opt=OFF smt_rear_opt=OFF
	local targets=()
	if (($# > 11)); then targets=("${@:12}"); else targets=("${AMS_BINS[@]}"); fi
	[[ "$hyper" == "ON" ]] && hyper_opt=ON
	[[ "$weighted" == "ON" ]] && weighted_opt=ON
	if [[ "$smt_symb" -gt 0 ]]; then smt_opt=ON; smt_rear_opt=ON; fi

	cmake -S "$ROOT" -B "$ROOT/build" \
		-DCMAKE_BUILD_TYPE=Release \
		-DTHREAD_NUM="$threads" -DUSE_HYPERTHREADING="$hyper_opt" \
		-DNUMA_MASK_VAL="$mask" -DUSE_WEIGHTED_INTERLEAVING="$weighted_opt" \
		-DUSE_INTERLEAVING=OFF -DCFG_PATH="$ROOT/$cfg" \
		-DAMS_HIGH_THREAD_NUM="$high_threads" \
		-DAMS_SIZE=524288 -DAMS_MIN_NUM=32768 \
		-DAMS_PRT="$prt" -DAMS_ACC_THRES="$acc" -DAMS_DNR=ON \
		-DAMS_SYMB_THRES=66.0 -DAMS_NUMC_THRES=1.50 \
		-DSMTACT="$smt_opt" -DSMTACT_REAR="$smt_rear_opt" \
		-DSMTACT_QUANTILE=0.99 -DSMTACT_SYMB_THRES="$smt_symb" \
		-DSMTACT_NUMC_THRES="$smt_numc" \
		-DDN_MIN=256 -DGUD_ALPHA=0.5 -DGUD_MIN=256 -DHYB_THRES=0.95 \
		-DHASH_CONSTANT=2654435761 -DRUN_NUM="$run_num" \
		-DUSE_HUGE=OFF -DMEMSV=OFF -DMEMSVB=OFF -DMEM_MON=OFF \
		-DIN_STATS=OFF -DIN_VERIFY=OFF -DIN_INSP=OFF -DIN_DEBUG=OFF -DIN_EXAMINE=OFF
	cmake --build "$ROOT/build" --parallel "$BUILD_JOBS" --target "${targets[@]}"
}

run_ams_dataset() {
	local dataset=$1 cpu_set=$2 log_prefix=$3
	local bins=("${AMS_BINS[@]}")
	if (($# > 3)); then bins=("$4"); fi
	local org="$ROOT/dataset/${dataset}.org.mtx"
	local trs="$ROOT/dataset/${dataset}.trs.mtx"
	if [[ ! -s "$org" || ! -s "$trs" ]]; then
		echo "Missing prepared MatrixMarket pair for '$dataset': $org and $trs" >&2
		return 1
	fi
	for bin in "${bins[@]}"; do
		local log="$LOG_DIR/${log_prefix}_${bin}.log"
		echo "[$(date '+%F %T')] $bin $dataset -> $log"
		"$NUMACTL" --physcpubind="$cpu_set" "$ROOT/bin/$bin" "$dataset" >> "$log" 2>&1
	done
}
