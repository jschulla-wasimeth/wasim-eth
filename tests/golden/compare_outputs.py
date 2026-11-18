#!/usr/bin/env python3
"""Compare two WaSiM output directories file by file.

Used as a golden test: run a setup with a reference binary and with a candidate
binary, then compare the two output directories.

- netCDF files (*.nc): every numeric variable is compared element-wise.
- all other files are read as text and split into whitespace-separated tokens:
  numeric tokens are compared with a tolerance, other tokens must match exactly.
  Tokens containing a path separator are ignored, so runs in different
  directories can be compared.
- log and diagnostics files are skipped.

A value pair passes if |a - b| <= atol + rtol * |b| (as in numpy.isclose).
Both missing and extra files are reported as failures.

usage: compare_outputs.py REFERENCE_DIR CANDIDATE_DIR [--rtol R] [--atol A] [--verbose]
exit code: 0 if all files match within tolerance, 1 otherwise
"""

import argparse
import math
import os
import sys

SKIP_SUFFIXES = (".log", ".xml", ".vrt", ".intern", "run_stats.txt", "control.wasim")


def is_number(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def compare_text(ref_path, cand_path, rtol, atol):
    """Return (n_values, max_abs_diff, max_rel_diff, first_problem or None)."""
    with open(ref_path, "r", encoding="latin-1") as f:
        ref = f.read().split()
    with open(cand_path, "r", encoding="latin-1") as f:
        cand = f.read().split()
    if len(ref) != len(cand):
        return 0, math.inf, math.inf, f"token count differs ({len(ref)} vs {len(cand)})"
    n = 0
    max_abs = 0.0
    max_rel = 0.0
    problem = None
    for i, (a, b) in enumerate(zip(ref, cand)):
        if a == b:
            if is_number(a):
                n += 1
            continue
        if "/" in a or "\\" in a:
            continue
        if is_number(a) and is_number(b):
            fa, fb = float(a), float(b)
            n += 1
            if math.isnan(fa) and math.isnan(fb):
                continue
            d = abs(fa - fb)
            rel = d / abs(fa) if fa != 0 else (math.inf if d > 0 else 0.0)
            max_abs = max(max_abs, d)
            max_rel = max(max_rel, rel)
            if d > atol + rtol * abs(fa) and problem is None:
                problem = f"token {i}: {a} vs {b}"
        elif problem is None:
            problem = f"token {i}: '{a}' vs '{b}'"
            max_abs = math.inf
    return n, max_abs, max_rel, problem


def compare_netcdf(ref_path, cand_path, rtol, atol):
    import numpy as np
    import netCDF4

    n = 0
    max_abs = 0.0
    max_rel = 0.0
    problem = None
    with netCDF4.Dataset(ref_path) as r, netCDF4.Dataset(cand_path) as c:
        if set(r.variables) != set(c.variables):
            return 0, math.inf, math.inf, "variable names differ"
        for name, var in r.variables.items():
            a = np.ma.filled(var[:], np.nan).astype("f8")
            b = np.ma.filled(c.variables[name][:], np.nan).astype("f8")
            if a.shape != b.shape:
                return n, math.inf, math.inf, f"{name}: shape {a.shape} vs {b.shape}"
            both_nan = np.isnan(a) & np.isnan(b)
            d = np.where(both_nan, 0.0, np.abs(a - b))
            if np.isnan(d).any():
                problem = problem or f"{name}: nodata pattern differs"
                d = np.nan_to_num(d, nan=np.inf)
            n += a.size
            if a.size:
                max_abs = max(max_abs, float(d.max()))
                with np.errstate(divide="ignore", invalid="ignore"):
                    rel = np.where(a != 0, d / np.abs(a), np.where(d > 0, np.inf, 0.0))
                max_rel = max(max_rel, float(np.nanmax(rel)))
                bad = d > atol + rtol * np.abs(np.nan_to_num(a))
                if bad.any() and problem is None:
                    idx = np.unravel_index(np.argmax(bad), a.shape)
                    problem = f"{name}{list(idx)}: {a[idx]} vs {b[idx]}"
    return n, max_abs, max_rel, problem


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference")
    ap.add_argument("candidate")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-9)
    ap.add_argument("--verbose", "-v", action="store_true", help="also list files that match")
    args = ap.parse_args()

    ref_files = {f for f in os.listdir(args.reference) if not f.endswith(SKIP_SUFFIXES)}
    cand_files = {f for f in os.listdir(args.candidate) if not f.endswith(SKIP_SUFFIXES)}
    failures = 0
    for f in sorted(ref_files - cand_files):
        print(f"MISSING  {f}")
        failures += 1
    for f in sorted(cand_files - ref_files):
        print(f"EXTRA    {f}")
        failures += 1

    total_values = 0
    worst = (0.0, "")
    for f in sorted(ref_files & cand_files):
        rp = os.path.join(args.reference, f)
        cp = os.path.join(args.candidate, f)
        if os.path.isdir(rp):
            continue
        try:
            if f.endswith(".nc"):
                n, dabs, drel, problem = compare_netcdf(rp, cp, args.rtol, args.atol)
            else:
                n, dabs, drel, problem = compare_text(rp, cp, args.rtol, args.atol)
        except Exception as e:  # unreadable file counts as a failure
            n, dabs, drel, problem = 0, math.inf, math.inf, f"error: {e}"
        total_values += n
        if dabs > worst[0]:
            worst = (dabs, f)
        if problem:
            failures += 1
            print(f"DIFF     {f}: {problem} (max abs {dabs:.3g}, max rel {drel:.3g})")
        elif args.verbose:
            print(f"ok       {f}: {n} values, max abs {dabs:.3g}")

    print(f"{len(ref_files & cand_files)} files compared, {total_values} values, "
          f"largest absolute difference {worst[0]:.3g} ({worst[1] or '-'}), {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
