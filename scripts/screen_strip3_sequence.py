#!/usr/bin/env python3
"""Bounded shared-memo prefixes; compare cache scheduling without stored answers.

F = full strip; A = end-middle deletion; D = both-end-middle deletion.
References are compared AFTER the process runs, never passed to the solver.
Completed jobs are flushed by STRIP3_TARGET_BATCH_PROGRESS; no timeout result
is interpreted as an SG value. Defaults: 2 s/process, 20 s/batch, one core,
512 MiB. Non-F orders stop at the available auxiliary reference width.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile
import time

from screen_strip3 import known_values, limits

ROOT = Path(__file__).resolve().parents[1]


def references():
    result = {("F", n): g for n, g in known_values().items()}
    for line in (ROOT / "data/strip3_language_gaps_32.out").read_text().splitlines():
        n, a, d = map(int, line.split())
        result["A", n] = a; result["D", n] = d
    ext = json.loads((ROOT / "data/strip3_a_ext.json").read_text())
    result.update({("A", int(n)): g for n, g in ext["values"].items()})
    return result


def jobs(order, maximum, known):
    labels, rows = [], []
    for n in range(1, maximum+1):
        if any((kind, n) not in known for kind in order):
            break
        full = (1 << (3*n))-1
        for kind in order:
            mask = full if kind == "F" else full & ~(1 << (3*n-2))
            if kind == "D": mask &= ~2
            labels.append([kind, n]); rows.append(f"{n} {mask:x}\n")
    return labels, "".join(rows)


def measure(binary, labels, input_text, known, env, timeout):
    start = time.monotonic()
    child = subprocess.Popen([str(binary), "masks"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True, env=env, start_new_session=True,
                             preexec_fn=lambda: limits(512, timeout))
    try:
        stdout, stderr = child.communicate(input_text, timeout=timeout)
        status = "complete" if child.returncode == 0 else "failed"
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGKILL); stdout, stderr = child.communicate()
        status = "timeout"
    try:
        values = list(map(int, stdout.split()))
    except ValueError:
        values = []; status = "invalid_output"
    if len(values) > len(labels): values = []; status = "invalid_output"
    completed = [dict(kind=kind, n=n, SG=g, expected=known[kind, n])
                 for (kind, n), g in zip(labels, values)]
    if any(row["SG"] != row["expected"] for row in completed): status = "wrong_value"
    elif status == "complete" and len(values) != len(labels): status = "incomplete_output"
    stages = []
    for line in stderr.splitlines():
        if line.startswith("batch "):
            stage = {key: float(value) if key == "time" else int(value)
                     for key, value in re.findall(r"(\w+)=([0-9.]+)", line)}
            stages.append(stage)
    return dict(status=status, returncode=child.returncode, wall_s=round(time.monotonic()-start, 6),
                requested_jobs=len(labels), completed=completed, stages=stages,
                next_job=labels[len(values)] if len(values) < len(labels) else None,
                stdout=stdout.strip(), stderr=stderr.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-source", type=Path, default=ROOT / "src/cfast/strip3_target.c")
    parser.add_argument("--orders", choices=("F", "ADF", "DAF"), nargs="+", default=["F", "ADF", "DAF"])
    parser.add_argument("--max-n", type=int, default=34)
    parser.add_argument("--env", action="append", default=[])
    parser.add_argument("--per-run", type=float, default=2)
    parser.add_argument("--total", type=float, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.max_n <= 34 or any(not math.isfinite(t) or t <= 0 for t in (args.per_run, args.total)):
        parser.error("known prefix through 34 and positive finite budgets required")
    extra = {}
    for setting in args.env:
        name, sep, value = setting.partition("=")
        if not sep or not re.fullmatch(r"STRIP3_[A-Z0-9_]+", name): parser.error("--env requires STRIP3_NAME=VALUE")
        extra[name] = value
    cpu = min(os.sched_getaffinity(0)); os.sched_setaffinity(0, {cpu})
    source_paths = set((args.candidate_source.resolve(), ROOT / "src/cfast/strip3_target.c",
                        ROOT / "src/cfast/strip3.c", ROOT / "src/cfast/strip3_comb_tables.h", Path(__file__).resolve()))
    report = dict(format="strip3-sequence-screen-v1", status="running", samples=[],
                  configuration=dict(orders=args.orders, max_n=args.max_n, extra_env=extra),
                  limits=dict(per_process_s=args.per_run, total_s=args.total, cpu=cpu, threads=1, address_space_mib=512),
                  provenance=dict(commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}))
    start = time.monotonic()
    known = references()
    def save():
        report["wall_s"] = round(time.monotonic()-start, 6)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n")
    env = {key: value for key, value in os.environ.items() if not key.startswith("STRIP3_")}
    env.update(STRIP3_TARGET_BITS="22", STRIP3_TARGET_ORDER="2", STRIP3_TARGET_SHARE_OPTIONS="1",
               STRIP3_TARGET_CACHE_FIRST="1", STRIP3_TARGET_FRONTIER="4", STRIP3_THREADS="1")
    env.update(extra); env["STRIP3_TARGET_BATCH_PROGRESS"] = "1"
    with tempfile.TemporaryDirectory(prefix="strip3-sequence-") as folder:
        binary = Path(folder) / "solver"
        try:
            subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                            str(args.candidate_source)], check=True, capture_output=True, timeout=args.total)
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
            report.update(status="build_failed", error=str(error)); save(); return
        for order in args.orders:
            labels, input_text = jobs(order, args.max_n, known)
            remaining = args.total-(time.monotonic()-start)
            if remaining <= 0: report["status"] = "total_budget"; save(); return
            row = measure(binary, labels, input_text, known, env, min(args.per_run, remaining))
            row["order"] = order; report["samples"].append(row); save()
            print(json.dumps({k: row[k] for k in ("order", "status", "wall_s", "next_job")}
                             | dict(completed_jobs=len(row["completed"]), last_stage=row["stages"][-1:])), flush=True)
            if row["status"] not in ("complete", "timeout"): report["status"] = "candidate_failed"; save(); return
            # These are alternative schedules, each stops at its first unknown
            # job; no failed schedule resumes with a larger time budget.
    report["status"] = "screened_requested_orders"; save()


if __name__ == "__main__":
    main()
