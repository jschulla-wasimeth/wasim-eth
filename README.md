# WaSiM, the Water balance Simulation Model

This repository contains the sources of the hydrologic model system WaSiM-ETH (Water balance Simulation Model, originally developed at ETH) forked from version 11.01.01. Earlier versions are available via https://www.wasim.ch/en/dialog/contact.htm

The model description (manual) is in [docs/wasim_2025_en.pdf](docs/wasim_2025_en.pdf), also available from https://www.wasim.ch/downloads/doku/wasim/wasim_2025_en.pdf.

Further documentation in [docs/](docs/):

- [sample_control_file.wasim](docs/sample_control_file.wasim): annotated sample control file with references to the manual
- [domain-size-limits.md](docs/domain-size-limits.md): limits on model domain size and memory
- [glossary.md](docs/glossary.md): German terms and abbreviations used in the code and file names

Known problems and planned work are listed in [TODO.md](TODO.md). Golden (regression) tests are described in [tests/golden/README.md](tests/golden/README.md).

## Dependencies

Required:

- a C++11 compiler (g++ 11 or newer on Linux, Apple clang on macOS)
- OpenMP (part of g++; on macOS provided by `libomp`)
- netCDF 4.x C library with headers (e.g. 4.6.1). The netCDF C++ interface is bundled in `src/netcdfcpp_Linux` and does not need to be installed separately.
- make

Optional:

- OpenMPI, for the parallel MPI builds (makefiles with `openmpi` / `mpi` in their name). Without it, use one of the non-MPI makefiles, which build an OpenMP-only binary.

Package names for common systems:

| System | Packages |
| --- | --- |
| RHEL / CentOS | `gcc-c++ make netcdf-devel` (EPEL), optional `openmpi-devel` |
| Ubuntu / Debian | `g++ make libnetcdf-dev`, optional `libopenmpi-dev openmpi-bin` |
| openSUSE | `gcc-c++ make netcdf-devel`, optional `openmpi4-devel` |
| macOS (Homebrew) | `xcode-select --install`, then `brew install netcdf llvm` (or `libomp` for Apple clang), optional `brew install open-mpi` |

## Compiling

`src/Makefile` builds on Linux and macOS:

    cd src
    make                  # OpenMP build      -> src/wasim
    make MPI=1            # MPI + OpenMP      -> src/wasim-mpi
    make NCXX=system      # use the system netCDF C++ library (ncxx4-config) instead of the bundled copy
    make help             # all options

The default build is portable. `make ARCH=-march=native` tunes it for the build machine, but this changes results slightly (fused multiply-add), so golden tests then need a looser tolerance. Objects go to `src/build/`, and `make clean` removes them along with the products of the older makefiles.

The per-distribution makefiles (`src/makefile_*`) are kept for now but superseded. If you still use them on Linux, note that they expect the netCDF C++ sources in `src/netcdfcpp` (`ln -s netcdfcpp_Linux netcdfcpp`) and list two source files that no longer exist (see TODO.md).

On macOS, install `netcdf llvm` (and `open-mpi` for MPI) with Homebrew. The Makefile then builds with Homebrew's LLVM clang, which has OpenMP built in, like clang on Linux; with MPI=1, Open MPI's `mpicxx` is pointed at that compiler. Without `llvm`, or with `make APPLE_CLANG=1`, it uses Apple clang and needs `brew install libomp`. Paths come from `brew --prefix` and can be overridden (`LLVM=...`, `LIBOMP=...`). `src/makefile_macos` (Apple clang only) is kept for compatibility.

## Running

    src/wasim path/to/control_file

WaSiM writes a translated copy of the control file (`<control file>.intern`), some routing tables to the current directory, and all other output to the `DefaultOutputDirectory` given in the control file. Use absolute paths for the output directory: with a relative one, statistics files of multi-layer land use end up in a doubled path (see TODO.md).

### Threads (OpenMP)

`OMP_NUM_THREADS` sets the number of threads. The parallel loops use the runtime schedule, so `OMP_SCHEDULE` applies to them. The default static schedule leaves threads idle when grid rows differ in work. `OMP_SCHEDULE=dynamic,4` made a 10-day winter run of the Rhone sample 27% faster on 4 threads, but a full year about twice as slow, so keep the default until this is understood (see TODO.md).

Results are bit-reproducible only with one thread. With several threads, a few cells differ between runs (see TODO.md).

### MPI

The MPI build splits the active cells of the domain into one region per process (rank), so a model can use several nodes. Rank 0 reads the input grids, gathers results for output and does the river routing. Each rank can additionally use OpenMP threads.

    OMP_NUM_THREADS=1 mpirun -n 4 src/wasim-mpi path/to/control_file      # srun on Slurm clusters such as Euler

- **When it helps:** for single runs that are too slow on one node. For ensembles (e.g. calibration), many independent single-node runs use the cores better: for the Doubs sample, 4 ranks gave a 2.2× speedup.
- **Memory:** rank 0 still needs the whole domain in memory; the other ranks need only their share plus about 110 MB.
- **Reproducibility:** at start-up each rank times a short benchmark, and faster ranks get more cells, so the decomposition can change between runs. A file `MpiSpeedFactor.dat` in the working directory, holding one number (the benchmark time, identical for all ranks), skips the benchmark and fixes the split.
- **Not available with MPI:** re-infiltration from rivers into groundwater (manual 2.1, p. 18).

### Debugging

Build with `make DEBUG=1` (`-O0 -g`, not stripped) and run under `gdb` or the VS Code C++ debugger with the control file as program argument. Use `OMP_NUM_THREADS=1` while debugging.
