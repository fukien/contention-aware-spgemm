# Taming Multi-Dimensional Skew in Sparse Matrix Multiplication with Contention-Aware Scheduling

This repository contains the source code for the paper *Taming Multi-Dimensional Skew in Sparse Matrix Multiplication with Contention-Aware Scheduling*.

[![arXiv preprint](https://img.shields.io/badge/arXiv-2610.36822-b31b1b.svg)](https://arxiv.org/abs/2610.36822)
[![Status](https://img.shields.io/badge/status-accepted-brightgreen)](https://arxiv.org/abs/2610.36822)

The preprint describes Adaptive Morsel Scheduling (AMS), including contention-aware scheduling and a hybrid accumulation method. It accompanies the paper accepted at VLDB 2027. See the [arXiv preprint](https://arxiv.org/abs/2610.36822) for the abstract and paper.

This codebase builds three SpGEMM implementations: `ams_hashspgemm`, `ams_hserscspgemm`, and `ams_hybaccspgemm`. It also builds the R-MAT conversion tools `rmat_ones` and `rmat_load`; PaRMAT can be built optionally.

---

## Prerequisites

The code targets Linux x86-64. The supplied machine configuration describes a two-socket, 64-core machine. For experiments using the original memory placement, the intended system has 4th-generation Intel Xeon Scalable processors, Linux kernel 6.9 or later, and CXL memory. The algorithms can also be built and run on NUMA DRAM systems; CXL is required only for CXL experiments.

Install the build dependencies on Debian or Ubuntu:

```bash
sudo apt update
sudo apt install build-essential cmake libnuma-dev libconfig-dev libpmem-dev
```

CMake 3.18 or later, a C11/C++17 compiler, OpenMP, pthreads, libnuma, libconfig, and libpmem are required. MKL, VTune, and PAPI are not used by this codebase.

## Machine configuration

The default CMake configuration uses `config/mc/mxc.cfg`, with 64 total cores and 32 cores per socket. Edit `config/mc/mxc.cfg` and `config/mc/mxc0.cfg` to match the core IDs and local/remote socket arrangement on your machine. Set `CORE_NUM`, `SOCKET_CORE_NUM`, and `CFG_PATH` when building if your topology differs.

For example, `mxc.cfg` treats socket 1 as local and socket 0 as remote. `mxc0.cfg` reverses that mapping. `CFG_PATH` defaults to `config/mc/mxc.cfg` and can be changed at configure time.

## Build

Build the three algorithms and two R-MAT utilities:

```bash
./build.sh
```

Add `-DBUILD_PARMAT=ON` to build the bundled R-MAT generator too. For a clean rebuild, run `./refresh.sh`; `./clean.sh` removes the generated `bin/` and `build/` directories.

The build scripts accept CMake options. `BUILD_JOBS` sets build parallelism (default 4). For example:

```bash
./build.sh -DTHREAD_NUM=32 -DAMS_HIGH_THREAD_NUM=4 -DRUN_NUM=2 \
  -DCFG_PATH=config/mc/mxc.cfg -DNUMA_MASK_VAL=1
```

`RUN_NUM` includes one warmup iteration and must be at least 2. `AMS_HIGH_THREAD_NUM` must be positive and smaller than `THREAD_NUM`. See `cmake/Options.cmake` for algorithm and memory-placement settings, including `SMTACT`, `SMTACT_REAR`, `AMS_DNR`, and AMS thresholds. NUMA masks use zero-based node bits; `LOCAL_NUMA_NODE_IDX` and `CXL_NUMA_NODE_IDX` use node ID plus one.

## Quick Start

Build PaRMAT and generate a symmetric R-MAT matrix. The converter expects a filename of the form `<prefix>_<scale>_<edge_factor>.txt`, and expects loop-free, duplicate-free, undirected input. For instance:

```bash
./build.sh -DBUILD_PARMAT=ON
mkdir -p dataset/rmat
./bin/PaRMAT -nVertices 262144 -nEdges 4194304 -threads 8 \
  -a 0.57 -b 0.19 -c 0.19 \
  -sorted -noEdgeToSelf -noDuplicateEdges -undirected \
  -output dataset/rmat/g500_18_16.txt
./bin/rmat_ones g500_18_16
./bin/rmat_load g500_18_16
```

`rmat_ones` writes `dataset/rmat/g500_18_16.org.mtx` and `dataset/rmat/g500_18_16.trs.mtx`. It sets matrix values to one and creates the transposed matrix. The SpGEMM executables read these prepared files; run them from the repository root:

```bash
./bin/ams_hybaccspgemm g500_18_16
```

Adjust CPU pinning to match your machine. The supplied utilities prepare R-MAT data only. For other datasets, provide the original and transposed MatrixMarket files in the path expected by the selected input mode. The former SuiteSparse and TAMU preparation utilities are not included.

## Repository layout

- `src/`: three SpGEMM entry points, R-MAT utilities, and shared code.
- `src/inc/`: common types, machine configuration, and memory allocation.
- `src/mm/`: shared sparse matrix operations and SpGEMM helpers.
- `src/ams/`: adaptive morsel scheduling and hybrid accumulation helpers.
- `config/mc/`: example two-socket CPU topology files.
- `PaRMAT/`: R-MAT graph generator.

## Source provenance

The R-MAT utilities, common infrastructure, configuration files, and PaRMAT came from CASMM. The three AMS implementations and their helpers came from UJMM-HMEM and were renamed to use the `ams_` prefix. This repository contains no experiment scripts or the original paper's full reproduction workflow.

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
