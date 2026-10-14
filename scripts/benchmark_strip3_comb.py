#!/usr/bin/env python3
"""Serial local ablation of certified comb reductions on full empty boards.

Every sample is a fresh process with one thread, an identical preallocated
memo budget, and the same binary. All configurations must agree on SG and
repeat deterministically in state counts. No server is used.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "data/strip3_comb_benchmark.json"
FLAGS = (7, 15, 23, 31)


def fields(line):
    return {key: int(value) for key, value in re.findall(r"(\w+)=(\d+)\b", line)}


def measure(binary, n, flags, bits, timeout, columns=None, label=None):
    began = time.perf_counter()
    command = [str(binary), "board", str(n)] if columns is None else [str(binary), "masks"]
    input_text = None if columns is None else f"{n} {sum(bits << (3*c) for c,bits in enumerate(columns)):x}\n"
    result = subprocess.run(command, input=input_text, text=True, capture_output=True, timeout=timeout,
                            env={**os.environ, "STRIP3_REDUCTIONS": str(flags), "STRIP3_THREADS": "1",
                                 "STRIP3_HT_BITS": str(bits), "STRIP3_HT_PREALLOC": "1",
                                 "STRIP3_HT_LOAD": "70", "STRIP3_PROGRESS": "0"})
    elapsed = time.perf_counter() - began
    if result.returncode:
        raise RuntimeError(f"incomplete n={n} flags={flags}: {result.stderr}")
    record = next(line for line in result.stderr.splitlines() if line.startswith("stats "))
    done = next(line for line in result.stderr.splitlines() if line.startswith("done "))
    output = list(map(int, result.stdout.split()))
    if (columns is None and (output[:2] != [3, n] or len(output) != 3)) or (columns is not None and len(output) != 1):
        raise RuntimeError("incomplete board output")
    return {"workload": label or f"full_3x{n}", "n": n, "reductions": flags, "SG": output[-1], "wall_s": round(elapsed, 6),
            "solver_s": float(re.search(r"done ([\d.]+)s", done)[1]),
            **fields(done), **fields(record)}


def summarize(samples):
    summary = []
    for workload, n in sorted({(s["workload"], s["n"]) for s in samples}):
        reference = None
        for flags in FLAGS:
            group = [s for s in samples if s["workload"] == workload and s["reductions"] == flags]
            stable = {key: group[0][key] for key in ("SG", "memo", "calls", "moves", "single_comb", "chain_comb", "chain_formula")}
            if not all(all(s[key] == value for key, value in stable.items()) for s in group):
                raise RuntimeError("state counts or SG changed between repetitions")
            row = {"workload": workload, "n": n, "reductions": flags, "repetitions": len(group), **stable,
                   "median_solver_s": statistics.median(s["solver_s"] for s in group),
                   "median_wall_s": statistics.median(s["wall_s"] for s in group),
                   "min_wall_s": min(s["wall_s"] for s in group), "max_wall_s": max(s["wall_s"] for s in group)}
            if reference is None:
                reference = row
            if row["SG"] != reference["SG"]:
                raise RuntimeError("SG differs between reduction settings")
            row["memo_saved"] = reference["memo"] - row["memo"]
            row["memo_saved_pct"] = 100 * row["memo_saved"] / reference["memo"]
            # The existing CLI reports packed entries; the separate wide-key
            # table (121..126 bits) is not included. Do not call that a total
            # state saving for the sparse 42-column family benchmarks.
            row["memo_scope"] = "all entries" if n <= 40 else "packed entries only; wide-key entries omitted"
            if n > 40:
                row["memo_saved"] = row["memo_saved_pct"] = None
            row["moves_saved_pct"] = 100 * (reference["moves"] - row["moves"]) / reference["moves"]
            row["wall_speedup"] = reference["median_wall_s"] / row["median_wall_s"]
            summary.append(row)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[18, 20, 22, 24])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--ht-bits", type=int, default=24)
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    args = parser.parse_args()
    if not args.sizes or any(n < 1 or n > 42 for n in args.sizes) or args.repetitions < 1:
        parser.error("positive repetitions and board sizes 1..42 required")
    cpu = min(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    if cpu is not None:
        os.sched_setaffinity(0, {cpu})
    provenance = {"platform": platform.platform(), "cpu_affinity": cpu,
                  "threads": 1, "ht_bits": args.ht_bits, "preallocated": True,
                  "compiler_flags": ["-O3", "-march=native", "-Wall", "-Wextra", "-pthread"],
                  "compiler": subprocess.check_output(["gcc", "--version"], text=True).splitlines()[0],
                  "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  "source_sha256": hashlib.sha256((ROOT / "src/cfast/strip3.c").read_bytes()).hexdigest(),
                  "tables_sha256": hashlib.sha256((ROOT / "src/cfast/strip3_comb_tables.h").read_bytes()).hexdigest()}
    provenance["memo_definition"] = "C CLI memo counts packed 120-bit keys; total state savings reported only for n<=40"
    samples = []
    with tempfile.TemporaryDirectory(prefix="strip3-comb-bench-") as folder:
        binary = Path(folder) / "strip3"
        subprocess.run(["gcc", *provenance["compiler_flags"], "-o", str(binary), str(ROOT / "src/cfast/strip3.c")], check=True)
        provenance["binary_sha256"] = hashlib.sha256(binary.read_bytes()).hexdigest()
        for n in args.sizes:
            for repeat in range(args.repetitions):
                # Rotate order to reduce a fixed first/last-position bias.
                offset = repeat % len(FLAGS)
                for flags in FLAGS[offset:] + FLAGS[:offset]:
                    sample = measure(binary, n, flags, args.ht_bits, args.timeout)
                    sample["repeat"] = repeat
                    samples.append(sample)
                    print(json.dumps(sample), flush=True)
            result = {"format": "node-kayles-comb-search-ablation-v1", "status": "partial",
                      "scope": "full cores K2..K5 and chain K2,K3 only; no damaged-core recognition",
                      "provenance": provenance, "samples": samples, "summary": summarize(samples)}
            args.output.write_text(json.dumps(result, indent=2) + "\n")
        arm = lambda length: tuple(4 if i % 2 == 0 else 6 for i in range(length))
        targets = {
            "single_K3_18_19": arm(18)[::-1] + (6,7,7,7,6) + arm(19),
            "chain_K2_K3_11_11_11": arm(11)[::-1] + (6,7,7,6) + arm(11) + (6,7,7,7,6) + arm(11),
        }
        for label, columns in targets.items():
            for repeat in range(args.repetitions):
                offset = repeat % len(FLAGS)
                for flags in FLAGS[offset:] + FLAGS[:offset]:
                    sample = measure(binary, len(columns), flags, args.ht_bits, args.timeout, columns, label)
                    sample["repeat"] = repeat
                    samples.append(sample)
                    print(json.dumps(sample), flush=True)
    result.update(status="complete", samples=samples, summary=summarize(samples))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "summary": result["summary"]}), flush=True)


if __name__ == "__main__":
    main()
