# Golden tests

Regression tests that run real model setups and compare every output file against reference outputs made with an earlier binary. They catch any change in results, intended or not, after refactoring, build changes or bug fixes.

The setups and case definitions are too large for this repository. They live in [wasim-forclim-evoland-experiments](https://github.com/mmyrte/wasim-forclim-evoland-experiments), under `golden/` (templates) and as git-lfs archives (setups):

| Case | Setup archive | Period | Active cells | Run time, 1 thread |
| --- | --- | --- | --- | --- |
| doubs | `sample-run-doubs.zip` | 2021-01-01 to 01-11, hourly | 6,022 | about 20 s |
| rhone | `sample-run-rhone.zip` | 2021-01-01 to 01-11, hourly | 21,584 | about 50 s |

Each case has a control file template (`golden/<case>/control.wasim.in`) with path placeholders; its header says how it differs from the original setup. The reference outputs are in `golden/<case>/reference.zip`; `golden/README.md` there says how they were made (single-threaded, from commit `a9b28ae`).

## Running

```sh
# in a checkout of wasim-forclim-evoland-experiments
unzip sample-run-rhone.zip -d cases/
python3 golden/rhone/merge_station_files.py cases/sample-run-rhone
unzip golden/rhone/reference.zip -d golden/rhone

# run the candidate binary
python3 tests/golden/run_case.py --binary src/wasim \
    --template golden/rhone/control.wasim.in \
    --case-dir cases/sample-run-rhone --out runs/rhone-candidate

# compare
python3 tests/golden/compare_outputs.py golden/rhone/reference runs/rhone-candidate
```

`compare_outputs.py` needs `numpy` and `netCDF4` (for the `.nc` outputs). It exits with 1 if any file is missing or differs beyond the tolerance (default `rtol=1e-6`, `atol=1e-9`).

## Reproducibility

- **Single-threaded runs are bit-identical**, which is why `run_case.py` defaults to `--threads 1` and the references are single-threaded.
- **With several OpenMP threads, results vary slightly between runs**: up to about 2e-6 absolute in reservoir balance files (`spv_*.dat`) for Rhone, from the order of parallel sums. Compare multi-threaded runs only with a looser `--atol 1e-5`.
- Results also depend on compiler and flags. The references were built with clang 18 at `-O3` without `-march=native`; a different compiler may need a looser tolerance, which is itself useful information.

## Updating references

Only replace the references when a change in results is intended and understood, and record which commit produced them and why they changed.
