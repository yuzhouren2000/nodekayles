#!/usr/bin/env python3
"""Bounded connection-matrix rank experiments for bitwise SG and P/N.

Different right behaviours need not be linearly independent. Test the stronger
linear compression proposal with GF(2) rank, keeping nonsingular-minor witnesses.
The 'guarded' family has two unplayed guard columns at each end; only the
outermost columns must be full. EVERY concatenation comes from a regular outer
history (unlike arbitrary induced-mask joins).
Finite rank data are not an infinite lower bound or a closed transfer algorithm.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import tempfile
import time

from interface_classes import all_codes
from strip3_history_language import outer_mask, words

ROOT = Path(__file__).resolve().parents[1]


def gf2_witness(rows):
    """Rank, independent original row indexes, and pivot column indexes.

    Every inserted row is reduced by previous pivots only, so the selected
    original rows and all pivot columns form a nonsingular square submatrix.
    """
    basis, selected, columns = {}, [], []
    for i, row in enumerate(rows):
        while row:
            pivot = row.bit_length()-1
            if pivot not in basis:
                basis[pivot] = row
                selected.append(i); columns.append(pivot)
                break
            row ^= basis[pivot]
    return {"rank": len(basis), "rows": selected, "columns": columns}


def catalogue(width, guarded=False):
    if not guarded:
        return [mask for mask, _ in all_codes(width)]
    if width < 4:
        raise ValueError("guarded widths must be at least 4")
    return sorted(outer_mask((0, 0)+word+(0, 0)) for word in words(width-4))


def row_bits(values, predicate):
    data = bytearray((len(values)+7)//8)
    for i, value in enumerate(values):
        if predicate(value):
            data[i//8] |= 1 << (i%8)
    return int.from_bytes(data, "little")


def ranks(matrix):
    maximum = max(map(max, matrix), default=0)
    selectors = [(f"SG_bit_{bit}", lambda g, bit=bit: (g >> bit) & 1)
                 for bit in range(max(1, maximum.bit_length()))]
    selectors += [("P", lambda g: g == 0)]
    result = {}
    packed = []
    width = len(matrix[0]) if matrix else 0
    for name, test in selectors:
        vectors = [row_bits(row, test) for row in matrix]
        result[name] = gf2_witness(vectors)
        if name != "P":
            packed.append(vectors)
    # A single shared linear row feature map for all SG bits must accommodate
    # the horizontal concatenation of the bit matrices, not just their max rank.
    vectors = [sum(packed[bit][i] << (bit*width) for bit in range(len(packed))) for i in range(len(matrix))]
    result["SG_shared_features"] = gf2_witness(vectors)
    return result


def resource_limits(timeout):
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (math.ceil(timeout)+1,)*2)


def measure(binary, objects, contexts, width, timeout, directory):
    jobs = directory / "jobs.txt"
    out = directory / "values.txt"
    err = directory / "stats.txt"
    with jobs.open("w") as stream:
        for left in objects:
            for right, size in contexts:
                stream.write(f"{width+size} {left | (right << (3*width)):x}\n")
    env = {key: value for key, value in os.environ.items() if not key.startswith("STRIP3_")}
    env.update(STRIP3_HT_BITS="22", STRIP3_THREADS="1", STRIP3_REDUCTIONS="7", STRIP3_PROGRESS="0")
    start = time.monotonic()
    with jobs.open() as stdin, out.open("w") as stdout, err.open("w") as stderr:
        child = subprocess.Popen([str(binary), "masks"], stdin=stdin, stdout=stdout, stderr=stderr,
                                 env=env, start_new_session=True, preexec_fn=lambda: resource_limits(timeout))
        try:
            child.wait(timeout=timeout)
            status = "complete" if child.returncode == 0 else "failed"
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL); child.wait()
            status = "timeout"
    info = {"status": status, "returncode": child.returncode,
            "wall_s": round(time.monotonic()-start, 6), "stderr": err.read_text().strip(),
            "queries": len(objects)*len(contexts)}
    if status != "complete":
        return info, None
    values = [int(value) for value in out.read_text().split()]
    if len(values) != info["queries"]:
        info["status"] = "incomplete_output"
        return info, None
    columns = len(contexts)
    return info, [bytearray(values[i*columns:(i+1)*columns]) for i in range(len(objects))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", type=int, nargs="+", default=[3, 4, 5, 6, 7, 8])
    parser.add_argument("--contexts", type=int, default=7, help="maximum width of reachable right contexts")
    parser.add_argument("--guarded", action="store_true")
    parser.add_argument("--per-run", type=float, default=2)
    parser.add_argument("--total", type=float, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (args.widths != sorted(set(args.widths)) or not args.widths or min(args.widths) < (4 if args.guarded else 1)
            or max(args.widths) > 12 or args.contexts < (4 if args.guarded else 1) or args.contexts > 12
            or args.per_run <= 0 or args.total <= 0 or not math.isfinite(args.per_run) or not math.isfinite(args.total)):
        parser.error("require increasing small widths and positive finite budgets")
    cpu = min(os.sched_getaffinity(0)); os.sched_setaffinity(0, {cpu})
    start = time.monotonic()
    source_paths = [ROOT / p for p in ("src/cfast/strip3.c", "src/cfast/strip3_comb_tables.h",
                                      "scripts/strip3_connection_rank.py", "scripts/strip3_history_language.py",
                                      "scripts/interface_classes.py")]
    report = {"format": "strip3-connection-rank-v1", "status": "running", "samples": [],
              "provenance": {"commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, cwd=ROOT).strip(),
                             "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}},
              "configuration": {"widths": args.widths, "contexts": args.contexts, "guarded": args.guarded},
              "limits": {"per_process_s": args.per_run, "total_s": args.total, "address_space_mib": 512,
                         "cpu": cpu, "threads": 1,
                         "total_budget_checks": "between width runs; solver deadlines are hard"},
              "interpretation": "Measured finite ranks only; rank witnesses require checking the actual SG entries."}

    def save():
        report["wall_s"] = round(time.monotonic()-start, 6)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n")

    context_sizes = range(4 if args.guarded else 1, args.contexts+1)
    contexts = [(mask, size) for size in context_sizes for mask in catalogue(size, args.guarded)]
    with tempfile.TemporaryDirectory(prefix="strip3-rank-") as folder:
        directory = Path(folder); binary = directory / "solver"
        try:
            subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                            str(ROOT / "src/cfast/strip3.c")], check=True, timeout=args.total, capture_output=True)
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
            report["status"] = "build_failed"; report["error"] = str(error); save(); return
        for width in args.widths:
            remaining = args.total-(time.monotonic()-start)
            if remaining <= 0:
                report["status"] = "total_budget"; save(); return
            objects = catalogue(width, args.guarded)
            info, matrix = measure(binary, objects, contexts, width, min(args.per_run, remaining), directory)
            info.update(width=width, objects=len(objects), contexts=len(contexts))
            if matrix is not None:
                info["ranks"] = ranks(matrix)
                info["distinct_behaviours"] = len(set(map(bytes, matrix)))
                info["maximum_SG"] = max(map(max, matrix))
                # Compact exact row bytes retain all entries for later checking.
                info["objects_hex"] = [f"{mask:x}" for mask in objects]
                info["contexts_hex"] = [[f"{mask:x}", size] for mask, size in contexts]
                info["matrix_hex"] = [bytes(row).hex() for row in matrix]
            report["samples"].append(info); save()
            print(json.dumps({k: v for k, v in info.items() if k not in ("matrix_hex", "objects_hex", "contexts_hex", "ranks")}
                             | {"ranks": {name: witness["rank"] for name,witness in info.get("ranks",{}).items()}}), flush=True)
            if info["status"] != "complete":
                report["status"] = "screen_stopped"; save(); return
    report["status"] = "completed_requested_widths"; save()


if __name__ == "__main__":
    main()
