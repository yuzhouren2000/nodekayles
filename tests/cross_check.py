"""Cross-check the C solver against the independent Python solver.

Both implementations were written separately from the same specification
(bitboard + connected-component Sprague-Grundy recursion + canonical form).
Agreement over all boards up to MMAX x NMAX is strong evidence against
implementation bugs in either.

Run:  python3 tests/cross_check.py [MMAX] [NMAX]      (default 6 8)
"""

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

from nodekayles import Solver  # noqa: E402

C_BIN = os.path.join(ROOT, "src", "cfast", "kayles")


def c_table(mmax, nmax):
    out = subprocess.run(
        [C_BIN, "table", str(mmax), str(nmax)],
        capture_output=True, text=True, check=True,
    ).stdout
    vals = {}
    for line in out.splitlines():
        m, n, g = line.split("\t")
        vals[(int(m), int(n))] = int(g)
    return vals


def main():
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    print(f"cross-checking all boards within {mmax} x {nmax}")
    t0 = time.time()
    cvals = c_table(mmax, nmax)
    print(f"  C solver: {len(cvals)} boards in {time.time()-t0:.1f}s")
    s = Solver()
    bad = 0
    t0 = time.time()
    for (m, n), cg in sorted(cvals.items()):
        pg = s.grid_value(m, n)
        if pg != cg:
            bad += 1
            print(f"  MISMATCH {m}x{n}: C={cg} Python={pg}")
    print(f"  Python solver: checked in {time.time()-t0:.1f}s, mismatches={bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
