# Model domain size: limits and overflow risks

Findings from reading the source of version 11.01.01, compared with the model description in [wasim_2025_en.pdf](wasim_2025_en.pdf) (manual version 11.00.10, June 2026). Nothing here was tested with a real large domain. The thresholds come from the data types in the code.

**Update 2026-10-11:** the netCDF buffer overflow, the zone code off-by-one errors and the `int` coordinates in the netCDF reader are fixed (wasim-eth commit history, see TODO.md). The sections below describe the original behaviour.

## Short answer

| Limit | Value | Where | Severity |
| --- | --- | --- | --- |
| Cells per grid, netCDF input/output | < 2,147,483,647 (≈ 46,341 × 46,341) | `int` products `ncols*nrows` in `modefuzr.cc`, `netcdfhlp.cc` | **fixed**: sizes now `size_t` (was: crash or memory corruption above it) |
| Cells per grid, Windows builds | ≈ 268 million (≈ 16,384 × 16,384) | `long` (32 bit on Windows) in `memory_handling::matrix()` | wrong allocation above it |
| Rows with MPI | 32,767 guaranteed by the MPI standard; Open MPI allows ~2^31 | row number used as MPI message tag | only matters for other MPI implementations |
| Zone codes (subbasin grid) | 0 … 32,767 (manual: up to 32,768) | `maxzonenwert`, `modconst.h` | **fixed** (was: code 32767 wrote out of bounds) |
| Number of distinct zones | 16,000 | `maxzonenzahl`, `modconst.h` | **fixed** (was: 16,000 zones wrote out of bounds) |
| Number of regions (regression interpolation) | 30 | `_max_regionen_zahl_` | checked |
| Grids held in memory | 1,000 | `gridmaxanzahl` | checked |
| Memory on MPI rank 0 | at least 1–2 full grids (8 bytes/cell each) during output | full-grid gather in `MPIGatherIntoNewGrid`, plus netCDF flat buffer | practical limit |

For ASCII and binary grids on Linux/macOS, nothing in the grid code overflows below roughly 2^31 rows or columns. In practice memory runs out first.

## Details

### How grids are stored

`memory_handling::matrix()` (`src/modefuzr.cc:1635`) allocates each grid as an array of row pointers, with the data in chunks of whole rows. Sizes are computed in `long`, which is 64 bit on Linux and macOS. A single row may be at most `_maxlinmem_` = 2,147,480,000 bytes, i.e. about 268 million columns. Cells are addressed as `grid[row][col]` with `int` row/column indices, so there are no flattened `row*ncols+col` indices in the model core.

Binary and ASCII grid files are read and written one row at a time (`readgridfromfile`, `printgridfile`), so file offsets do not limit grid size. Single-precision binary grids store the header (ncols, nrows, corner coordinates, cellsize) as `float`. That is exact for integers up to 16,777,216, which is enough for row/column counts and for metric corner coordinates. Corner coordinates are rounded down to whole metres when grids are compared anyway (`modefuzr.cc:2166ff`).

### netCDF I/O: `int` overflow (the main risk)

Whenever netCDF grids are written or read, the whole grid is copied into a flat buffer:

```cpp
double *gridbuf = new double[(ncols)*(nrows)];   // int * int
int count = 0;
```

(`src/modefuzr.cc:4063–4104`, `4358–4394`, `6014–6056`; `src/netcdfhlp.cc:880–960` also loops `for (int ii = 0; ii < nc_cols*nc_rows; ...)`).

`ncols*nrows` is computed in 32-bit `int`. Signed overflow is undefined behaviour; in practice it wraps:

- 50,000 × 50,000 = 2.5·10⁹ → wraps to a negative number, so `new[]` throws `std::bad_array_new_length` and the model aborts.
- 70,000 × 70,000 = 4.9·10⁹ → wraps to 605,032,704, so the buffer is too small and the copy loop writes far past its end (heap corruption, wrong results or a crash much later).

Grids of this size, for reference: 46,341 cells per side is about 1,160 km at 25 m, 463 km at 10 m, or 46 km at 1 m. Models without netCDF grids are not affected.

Fix: use `size_t` for the buffer size and the loop counters.

### Windows: 32-bit `long`

On 64-bit Windows `long` is 32 bit. In `matrix()`, `size = nrow*ncol*sizeof(double)` is stored in a `long`, so grids larger than 2^31 bytes (≈ 268 million cells) get a wrong size. Linux and macOS (64-bit `long`) are not affected.

### MPI

- The domain is split by cells (`defltuzr.cc:8540ff`, counted in `long`). Each rank holds its own rows, and halos are exchanged row by row.
- Row segments are sent with `MPI_Send(..., cnt, MPI_DOUBLE, ..., tag = row, ...)`. `cnt` is at most `ncols` (times `levels + 1` for soil stacks), well below the `int` count limit for realistic grids.
- The message tag is the row number. The MPI standard only guarantees tags up to 32,767 (`MPI_TAG_UB`). Open MPI supports up to 2^31−1, but other implementations may have lower limits. With more than 32,767 rows on such an implementation, sends fail with `MPI_ERR_TAG`.
- For output, rank 0 gathers each grid into a full-size grid (`MPIGatherIntoNewGrid`, `modefuzr.cc:3070`); with netCDF it also builds the flat buffer above. Rank 0 therefore needs memory for whole grids even when the computation is spread out. Budget at least 16 bytes per cell on rank 0 on top of its own share.

### Zone and region codes (off by one)

`src/defltuzr.cc:1583–1616`:

```cpp
zzarray = calloc(maxzonenwert, ...);          // indices 0 … 32766
if ((zonenid<0) || (zonenid>maxzonenwert))    // lets 32767 through
...
zzone = calloc(maxzonenzahl, ...);            // indices 0 … 15999
if (zonmax > maxzonenzahl)                    // lets 16000 through
```

Code 32767, or exactly 16,000 distinct zones, writes one element past the end of the array. The manual ([wasim_2025_en.pdf](wasim_2025_en.pdf), section 7.3 "Error codes of WaSiM") documents error 3 as "zone codes or region codes too small (<1) or too large (>32768)" and error 4 as "too many zones (>16000)", and gives pour points the range 1…32767. So values the manual allows can trigger the bug, and code 32768 is in fact rejected. The same pattern is used for regions (`defltuzr.cc:3783ff`). Until this is fixed, keep zone codes ≤ 32,766 and use fewer than 16,000 zones. Also, `zonenid = (int)zn[i][j] + 0.5` truncates rather than rounds, so codes must be stored as exact integers in the grid.

### netCDF coordinates

The netCDF reader stores the `x` and `y` coordinate variables in `int` arrays (`netcdfhlp.cc:819`, `864`) and derives `xll`, `yll` and, if there is no `default_cellsize` attribute, the cellsize from them. Grids with fractional coordinates are therefore not supported. That includes geographic coordinates in degrees (cellsize becomes 0) and cells smaller than 1 m. `add_offset` for packed variables is also read into an `int`.

### Other

- Dates are compared as `year*535680 + month*44640 + …` minutes in `long`. This overflows a 32-bit `long` (Windows) only after year ~4000, so it does not matter.
- GIF input grids are limited to 65,535 × 65,535 pixels by the GIF format itself.
