#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
THREAD_NUM=${THREAD_NUM:-8}
EDGE_FACTORS=(2 4 8 16 32)

cd "$ROOT"
mkdir -p "$ROOT/dataset"

# Rebuild the project and PaRMAT before preparing the datasets.
"$ROOT/refresh.sh" -DBUILD_PARMAT=ON

for family in er g500; do
	if [[ "$family" == "er" ]]; then
		SCALES=(19 20 21)
		a=0.25
		b=0.25
		c=0.25
	else
		SCALES=(17 18 19)
		a=0.57
		b=0.19
		c=0.19
	fi

	for scale in "${SCALES[@]}"; do
		num_vertices=$((2 ** scale))
		for edge_factor in "${EDGE_FACTORS[@]}"; do
			num_edges=$((edge_factor * num_vertices))
			dataset="${family}_${scale}_${edge_factor}"
			edgelist="$ROOT/dataset/${dataset}.txt"
			matrix_org="$ROOT/dataset/${dataset}.org.mtx"
			matrix_trs="$ROOT/dataset/${dataset}.trs.mtx"

			# Remove prior outputs so this run regenerates each dataset from scratch.
			rm -f -- "$edgelist" "$matrix_org" "$matrix_trs"

			"$ROOT/bin/PaRMAT" \
				-nVertices "$num_vertices" \
				-nEdges "$num_edges" \
				-output "$edgelist" \
				-a "$a" -b "$b" -c "$c" \
				-threads "$THREAD_NUM" \
				-sorted -noEdgeToSelf -noDuplicateEdges -undirected \
				-memUsage 0.4

			"$ROOT/bin/rmat_ones" "$dataset"
		done
	done
done

printf 'Prepared ER and G500 R-MAT datasets in %s/dataset\n' "$ROOT"
