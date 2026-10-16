#!/usr/bin/env python3
"""Fast local screening of a board-N CLI against established exact SG values.

Hard wall/CPU/address-space limits; one CPU; fresh processes; stop at the first
candidate failure. This is a screening experiment, not a scaling proof or a
repeated timing benchmark. No server, SSH, or long background job is involved.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import signal
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def limits(memory_mib, seconds):
    resource.setrlimit(resource.RLIMIT_AS, (memory_mib * 1024**2,) * 2)
    resource.setrlimit(resource.RLIMIT_CPU, (math.ceil(seconds) + 1,) * 2)


def measure(binary, n, timeout, memory_mib, env):
    start = time.monotonic()
    process = subprocess.Popen([str(binary), "board", str(n)], text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               env={**{key: value for key, value in os.environ.items()
                                       if not key.startswith("STRIP3_")}, **env}, start_new_session=True,
                               preexec_fn=lambda: limits(memory_mib, timeout))
    try:
        stdout, stderr = process.communicate(timeout=timeout)
        status = "complete" if process.returncode == 0 else "failed"
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        status = "timeout"
    row = {"n": n, "status": status, "wall_s": round(time.monotonic()-start, 6),
           "returncode": process.returncode, "stdout": stdout.strip(), "stderr": stderr.strip()}
    row["counts"] = {key: int(value) for key, value in re.findall(r"\b(\w+)=(\d+)\b(?!\.)", stderr)}
    obligations = []
    for line in stderr.splitlines():
        if line.startswith("obligation "):
            found = dict(re.findall(r"(\w+)=(-?[0-9a-f]+)", line))
            obligations.append({key: value if key == "mask" else int(value) for key, value in found.items()})
    if obligations:
        row["obligations"] = obligations
    fields = stdout.split()
    if status == "complete":
        if len(fields) != 3 or fields[:2] != ["3", str(n)]:
            row["status"] = "invalid_output"
        else:
            try:
                row["SG"] = int(fields[2])
            except ValueError:
                row["status"] = "invalid_output"
    # Incomplete rows are never treated as values, even if stdout has a prefix.
    return row


def known_values():
    values = {int(n): int(g) for _, n, g in
              (line.split() for line in (ROOT / "data/raw/strip3_3xn.tsv").read_text().splitlines())}
    n, g = map(int, (ROOT / "data/strip3_language_34.out").read_text().split())
    values[n] = g
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-source", type=Path, default=ROOT / "src/cfast/strip3_target.c")
    parser.add_argument("--order", type=int, choices=(0, 1, 2, 3, 4), default=2)
    parser.add_argument("--low-first", type=int, choices=(0, 1), default=0)
    parser.add_argument("--skip-reference", action="store_true")
    parser.add_argument("--env", action="append", default=[], metavar="NAME=VALUE",
                        help="record an additional STRIP3_* experiment setting")
    parser.add_argument("--sizes", type=int, nargs="+", default=[12, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34])
    parser.add_argument("--per-run", type=float, default=2)
    parser.add_argument("--total", type=float, default=20)
    parser.add_argument("--memory-mib", type=int, default=512)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    extra_env = {}
    for setting in args.env:
        name, separator, value = setting.partition("=")
        if not separator or not re.fullmatch(r"STRIP3_[A-Z0-9_]+", name):
            parser.error("--env requires STRIP3_NAME=VALUE")
        extra_env[name] = value
    known = known_values()
    if (not args.sizes or args.sizes != sorted(set(args.sizes)) or any(n not in known for n in args.sizes)
            or not math.isfinite(args.per_run) or not math.isfinite(args.total)
            or args.per_run <= 0 or args.total <= 0 or args.memory_mib < 128):
        parser.error("require increasing known sizes, positive finite deadlines, and at least 128 MiB")
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    start = time.monotonic()
    result = {"format": "strip3-target100-screen-v1", "status": "running", "samples": [],
              "configuration": {"order": args.order, "low_first": args.low_first,
                                "sizes": args.sizes, "skip_reference": args.skip_reference,
                                "extra_env": extra_env},
              "limits": {"per_process_wall_s": args.per_run, "total_wall_s": args.total,
                         "address_space_mib": args.memory_mib, "cpu": cpu, "threads": 1},
              "provenance": {"commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                             "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
                                               (args.candidate_source.resolve(), ROOT / "src/cfast/strip3.c",
                                                ROOT / "src/cfast/strip3_target.c",
                                                ROOT / "src/cfast/strip3_comb_tables.h")}},
              "interpretation": "One sample per size; finite measurements do not establish asymptotic complexity.",
              "timing_caveat": "Other local work may run concurrently; compare operation counts before timing ratios."}

    def save():
        result["wall_s"] = round(time.monotonic()-start, 6)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")

    reference_active = not args.skip_reference
    with tempfile.TemporaryDirectory(prefix="strip3-screen-") as folder:
        binaries = {}
        sources = [("candidate", args.candidate_source)]
        if not args.skip_reference:
            sources.append(("reference", ROOT / "src/cfast/strip3.c"))
        for name, source in sources:
            binary = Path(folder) / name
            remaining = args.total - (time.monotonic()-start)
            if remaining <= 0:
                result["status"] = "total_budget"; save(); return
            try:
                subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                                str(source)], check=True, capture_output=True, text=True, timeout=remaining)
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
                result["status"] = "build_failed"; result["error"] = str(error); save(); return
            binaries[name] = binary
        for n in args.sizes:
            for name in ("candidate", "reference"):
                if name == "reference" and not reference_active:
                    continue  # do not repeatedly spend the timeout on the old algorithm
                remaining = args.total - (time.monotonic()-start)
                if remaining <= 0:
                    result["status"] = "total_budget"; save(); return
                env = {"STRIP3_TARGET_BITS": "22", "STRIP3_TARGET_ORDER": str(args.order),
                       "STRIP3_TARGET_LOW_FIRST": str(args.low_first),
                       "STRIP3_THREADS": "1", "STRIP3_HT_BITS": "24", "STRIP3_HT_PREALLOC": "1",
                       "STRIP3_REDUCTIONS": "7", "STRIP3_PROGRESS": "0"}
                env.update(extra_env)
                row = measure(binaries[name], n, min(args.per_run, remaining), args.memory_mib, env)
                row["solver"] = name
                if row["status"] == "complete" and row["SG"] != known[n]:
                    row["status"] = "wrong_value"
                result["samples"].append(row)
                print(json.dumps(row), flush=True)
                save()
                if name == "candidate" and row["status"] != "complete":
                    result["status"] = "candidate_stopped"; save(); return
                if row["status"] == "wrong_value":
                    result["status"] = "reference_disagrees"; save(); return
                # Reference timeout is allowed: a faster candidate can continue.
                if name == "reference" and row["status"] != "complete":
                    reference_active = False
    result["status"] = "completed_requested_sizes"
    save()


if __name__ == "__main__":
    main()
