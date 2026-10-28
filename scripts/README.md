# Scripts

These scripts prepare datasets and run the three AMS SpGEMM algorithms. Run them from the repository root.

## R-MAT

```bash
./scripts/rmat_preparation.sh
./scripts/ams_rmat_run.sh
```

The preparation script generates the ER and G500 matrices under `dataset/`. The run script builds and runs `ams_hashspgemm`, `ams_hserscspgemm`, and `ams_hybaccspgemm` on those matrices. It builds each algorithm separately because the R-MAT experiment uses different AMS and SMTACT tuning for HASH versus HSERSC/HYB. The default `RUN_NUM` is 10, including the warmup. Logs go to `logs/ams-rmat/`.

## Join-Aggregate

```bash
./scripts/join-aggregate_preparation.sh
./scripts/ams_join-aggregate_run.sh
```

The preparation script downloads the source data and creates the CLH0 (`chbc_s_uni`), BKGN (`bkgn`), and IMDb (`imdb_movie`) matrix pairs under `dataset/`. The run script builds and runs all three AMS algorithms on them using a shared local-node profile without weighted interleaving. The default `RUN_NUM` is 10, including the warmup. Logs go to `logs/ams-join-aggregate/`.

## Shared helper

`ams_run_common.sh` is sourced by the two run scripts. It contains their shared CMake build and dataset execution functions; you do not run it directly.
