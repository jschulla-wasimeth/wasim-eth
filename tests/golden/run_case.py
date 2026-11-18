#!/usr/bin/env python3
"""Run one golden-test case with a given WaSiM binary.

The case is a control file template in which absolute paths are replaced by
placeholders, so the same template works wherever the data is unpacked:

  @CASE_DIR@    the unpacked model setup (Input/, StateIni_*/ ...)
  @GOLDEN_DIR@  the directory holding the template and any extra input files
  @OUT_DIR@     the output directory of this run (created; must not exist yet)

The run uses @OUT_DIR@ as working directory, because WaSiM writes some files
(e.g. routing tables) to the current directory. Peak memory and run time are
printed at the end and written to @OUT_DIR@/run_stats.txt.

usage: run_case.py --binary BIN --template TEMPLATE --case-dir DIR --out DIR [--threads N] [--mpi N]
"""

import argparse
import os
import resource
import subprocess
import sys
import time


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--template", required=True)
    ap.add_argument("--case-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1,
                    help="OMP_NUM_THREADS (default 1: only single-threaded runs are bit-reproducible)")
    ap.add_argument("--mpi", type=int, default=0, help="run with mpirun -n N (MPI build only)")
    args = ap.parse_args()

    binary = os.path.abspath(args.binary)
    template = os.path.abspath(args.template)
    case_dir = os.path.abspath(args.case_dir)
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=False)

    with open(template, encoding="latin-1") as f:
        text = f.read()
    text = (text.replace("@CASE_DIR@", case_dir)
                .replace("@GOLDEN_DIR@", os.path.dirname(template))
                .replace("@OUT_DIR@", out_dir))
    control = os.path.join(out_dir, "control.wasim")
    with open(control, "w", encoding="latin-1") as f:
        f.write(text)

    cmd = [binary, control]
    if args.mpi:
        cmd = ["mpirun", "-n", str(args.mpi)] + cmd
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads))
    start = time.time()
    with open(os.path.join(out_dir, "run.log"), "w") as log:
        rc = subprocess.run(cmd, cwd=out_dir, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    # ru_maxrss is in KiB on Linux and in bytes on macOS
    rss_mb = usage.ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024)
    stats = (f"exit {rc}  wall {time.time() - start:.1f} s  cpu {usage.ru_utime + usage.ru_stime:.1f} s  "
             f"peak memory {rss_mb:.0f} MB  threads {args.threads}  mpi {args.mpi}")
    with open(os.path.join(out_dir, "run_stats.txt"), "w") as f:
        f.write(stats + "\n")
    print(stats)
    return rc


if __name__ == "__main__":
    sys.exit(main())
