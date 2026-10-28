# Taming Multi-Dimensional Skew in Sparse Matrix Multiplication with Contention-Aware Scheduling

This repository contains the source code for the paper *Taming Multi-Dimensional Skew in Sparse Matrix Multiplication with Contention-Aware Scheduling*.

[![arXiv preprint](https://img.shields.io/badge/arXiv-2610.36822-b31b1b.svg)](https://arxiv.org/abs/2610.36822)
[![Status](https://img.shields.io/badge/status-accepted-brightgreen)](https://arxiv.org/abs/2610.36822)

The preprint describes Adaptive Morsel Scheduling (AMS), including contention-aware scheduling and a hybrid accumulation method. It accompanies the paper accepted at VLDB 2027. See the [arXiv preprint](https://arxiv.org/abs/2610.36822) for the abstract and paper.

This codebase builds three SpGEMM implementations: `ams_hashspgemm`, `ams_hserscspgemm`, and `ams_hybaccspgemm`. It also builds the R-MAT converter `rmat_ones`; PaRMAT can be built optionally.

---

## Prerequisites

The code targets Linux x86-64. The supplied machine configuration is an example for a two-socket, 64-core system. The included runners use local NUMA placement and can run on NUMA DRAM systems. CXL memory is needed only if you configure a CXL placement experiment; the reference setup uses 4th-generation Intel Xeon Scalable processors and Linux kernel 6.9 or later.

Install the build dependencies on Debian or Ubuntu:

```bash
sudo apt update
sudo apt install build-essential cmake libnuma-dev libconfig-dev libpmem-dev numactl curl python3
```

CMake 3.18 or later, a C11/C++17 compiler, OpenMP, pthreads, libnuma, libconfig, and libpmem are required. MKL, VTune, and PAPI are not used by this codebase.

## Machine configuration

The default CMake configuration uses `config/mc/mxc0.cfg`, an example topology with 64 core IDs and 32 cores per socket. Treat it as a template: edit the core IDs and socket arrangement to match your system. Set `CORE_NUM`, `SOCKET_CORE_NUM`, and `CFG_PATH` when building if your topology differs.

`mxc0.cfg` maps cores 0–31 as local and 32–63 as remote. `CFG_PATH` defaults to `config/mc/mxc0.cfg` and can be changed at configure time.

The default `NUMA_MASK_VAL` is `1`, which selects NUMA node 0. The `mxc0.cfg` file maps core IDs; it does not set the NUMA mask. Workload scripts can pass their own memory-placement settings.

## Build

Build the three algorithms and R-MAT converter:

```bash
./build.sh
```

Add `-DBUILD_PARMAT=ON` to build the bundled R-MAT generator too. For a clean rebuild, run `./refresh.sh`; `./clean.sh` removes the generated `bin/` and `build/` directories.

The build scripts accept CMake options. `BUILD_JOBS` sets build parallelism (default 4). For example:

```bash
./build.sh -DTHREAD_NUM=32 -DAMS_HIGH_THREAD_NUM=4 -DRUN_NUM=10 \
  -DCFG_PATH=config/mc/mxc0.cfg -DNUMA_MASK_VAL=1
```

`RUN_NUM` includes one warmup iteration and must be at least 2. `AMS_HIGH_THREAD_NUM` must be positive and smaller than `THREAD_NUM`. See `cmake/Options.cmake` for algorithm and memory-placement settings, including `SMTACT`, `SMTACT_REAR`, `AMS_DNR`, and AMS thresholds. `NUMA_MASK_VAL` uses node bits, while `LOCAL_NUMA_NODE_IDX` and `CXL_NUMA_NODE_IDX` use node ID plus one.

## Quick Start

Build PaRMAT and generate a symmetric R-MAT matrix. The converter expects a filename of the form `<prefix>_<scale>_<edge_factor>.txt`, and expects loop-free, duplicate-free, undirected input. For instance:

```bash
./build.sh -DBUILD_PARMAT=ON
mkdir -p dataset
./bin/PaRMAT -nVertices 262144 -nEdges 4194304 -threads 8 \
  -a 0.57 -b 0.19 -c 0.19 \
  -sorted -noEdgeToSelf -noDuplicateEdges -undirected \
  -output dataset/g500_18_16.txt
./bin/rmat_ones g500_18_16
```

`rmat_ones` writes `dataset/g500_18_16.org.mtx` and `dataset/g500_18_16.trs.mtx`. It sets matrix values to one and creates the transposed matrix. The SpGEMM executables read these prepared files; run them from the repository root:

```bash
./bin/ams_hybaccspgemm g500_18_16
```

Adjust CPU pinning to match your machine. The included scripts prepare R-MAT and Join-Aggregate datasets. For other datasets, place the original and transposed MatrixMarket files in `dataset/` as `<name>.org.mtx` and `<name>.trs.mtx`.

## Dataset Preparation

### R-MAT datasets

To generate the ER and G500 R-MAT datasets from scratch, run:

```bash
./scripts/rmat_preparation.sh
```

This rebuilds PaRMAT and generates ER scales 19–21 and G500 scales 17–19, each at edge factors 2, 4, 8, 16, and 32, using 8 generator threads and the configured ER/G500 parameters. Set `THREAD_NUM` to change the generator thread count. The sweep contains 30 matrices and can require substantial time, memory, and disk space.

### Join-Aggregate datasets

To download and prepare the CLH0, BKGN, and IMDB Join-Aggregate matrices, run:

```bash
./scripts/join-aggregate_preparation.sh
```

The script downloads ClickBench `hits.tsv.gz`, GroupLens Book Genome, and the IMDb `title.basics` and `title.principals` tables. It builds the ClickBench user–URL relation and samples 850,000 users uniformly with seed 42 for `chbc_s_uni` (CLH0), converts Book Genome ratings to `bkgn` (BKGN), and builds the movie-title–credited-person incidence matrix `imdb_movie` (IMDB). The downloaded archives and the resulting general MatrixMarket files are kept directly under `dataset/`, which is the input directory expected by the current runners. Temporary processing files are removed after preparation. This preparation needs substantial free disk space and time; IMDb's datasets are provided for personal and non-commercial use.

## Running the AMS experiments

After preparing the matching datasets, run the workload scripts from the repository root:

```bash
./scripts/ams_rmat_run.sh
./scripts/ams_join-aggregate_run.sh
```

Both scripts build and run `ams_hashspgemm`, `ams_hserscspgemm`, and `ams_hybaccspgemm`, and write fresh logs under `logs/ams-rmat/` or `logs/ams-join-aggregate/`. The R-MAT runner covers the 30 ER/G500 matrices above, uses 64 threads, and runs ten iterations including warmup. The Join-Aggregate runner uses local-node placement without weighted interleaving and also defaults to 64 threads and ten iterations. Set `THREAD_NUM`, `RUN_NUM`, `AMS_HIGH_THREAD_NUM`, `BUILD_JOBS`, or `LOG_DIR` in the environment to override these defaults. See [`scripts/README.md`](scripts/README.md) for script usage. Before running on a different system, adjust the CPU pinning ranges and the example topology in `config/mc/mxc0.cfg` to match its processor and memory layout.

## Repository layout

- `src/`: three SpGEMM entry points, R-MAT utilities, and shared code.
- `src/inc/`: common types, machine configuration, and memory allocation.
- `src/mm/`: shared sparse matrix operations and SpGEMM helpers.
- `src/ams/`: adaptive morsel scheduling and hybrid accumulation helpers.
- `scripts/`: dataset preparation scripts plus AMS runners for R-MAT and Join-Aggregate workloads.
- `config/`: machine configuration, including the default two-socket topology in `mc/mxc0.cfg`.
- `PaRMAT/`: R-MAT graph generator.

## Citation

If this code or the AMS method is useful in your research, please cite the arXiv preprint:

```bibtex
@misc{huang2026taming,
  title         = {{Taming Multi-Dimensional Skew in Sparse Matrix Multiplication with Contention-Aware Scheduling}},
  author        = {Wentao Huang and Qiming Huang and Yunhong Ji and Mian Lu and Yuqiang Chen and Bingsheng He and Kian-Lee Tan},
  year          = {2026},
  eprint        = {2610.36822},
  archivePrefix = {arXiv},
  primaryClass  = {cs.DC},
  url           = {https://arxiv.org/abs/2610.36822}
}
```

Preprint: [arXiv:2610.36822](https://arxiv.org/abs/2610.36822).
