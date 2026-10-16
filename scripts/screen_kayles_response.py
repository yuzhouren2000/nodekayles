#!/usr/bin/env python3
"""Bounded local screen of semantic response reduction in column transfer.

The worker never reads known SG data. The parent validates only completed
rows after termination; budget failures have no SG for the unfinished width.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def worker(mode, max_n, pair_limit):
    from kayles_boundary_transfer import BoundaryTransfer, ThreeRowTransfer
    from kayles_response_equivalence import ResponseTransfer
    algebra = ResponseTransfer(pair_limit=pair_limit) if mode == "response" else BoundaryTransfer()
    model = ThreeRowTransfer(algebra)
    state = (0, 0)
    started = time.monotonic()
    for n in range(1, max_n+1):
        try:
            state = model.append(state, 0 if n == 1 else 7, 7)
        except RuntimeError as error:
            print(json.dumps(dict(n=n, status="incomplete", reason=str(error), counts=algebra.report())), flush=True)
            return 3
        g = algebra.value(state)
        counts = algebra.report() if mode == "response" else dict(forms=len(algebra.nodes)-1,
                                                                 value_evaluations=algebra._value.cache_info().misses)
        print(json.dumps(dict(n=n, SG=g, status="complete", counts=counts,
                              elapsed_s=round(time.monotonic()-started, 6))), flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=("strict", "response"))
    parser.add_argument("--max-n", type=int, default=100)
    parser.add_argument("--reference-max-n", type=int, default=10)
    parser.add_argument("--pair-limit", type=int, default=200000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args()
    if args.worker:
        return worker(args.worker, args.max_n, args.pair_limit)
    if args.output is None or not 1 <= args.reference_max_n <= args.max_n <= 100 or args.pair_limit < 1:
        parser.error("require output, 1 <= reference-max-n <= max-n <= 100 and positive pair limit")
    from screen_strip3 import atomic_write_json, known_values
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    started = time.monotonic()
    sources = [Path(__file__).resolve(), ROOT / "scripts/kayles_boundary_algebra.py",
               ROOT / "scripts/kayles_boundary_transfer.py", ROOT / "scripts/kayles_response_equivalence.py",
               ROOT / "scripts/screen_strip3.py", ROOT / "data/raw/strip3_3xn.tsv",
               ROOT / "data/strip3_language_34.out"]
    fingerprint = dict(sources={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources}, python=sys.version, executable=sys.executable,
                       max_n=args.max_n, reference_max_n=args.reference_max_n, pair_limit=args.pair_limit,
                       per_process_wall_s=2, total_wall_s=20, address_space_mib=512, cpu=cpu)
    known = known_values()
    if args.output.exists() and not args.fresh:
        old = json.loads(args.output.read_text())
        if old.get("fingerprint") == fingerprint and old.get("status") in ("complete", "candidate_stopped"):
            for run in old["runs"]:
                for row in run["rows"]:
                    if row["status"] == "complete":
                        if row["n"] in known:
                            assert row["SG"] == known[row["n"]]
                    else:
                        assert "SG" not in row
            print(json.dumps(dict(evidence_reused=True, solver_processes=0, status=old["status"],
                                  note="reuse saved bounded experiment; no new SG or solver speed claim")))
            return 0
    report = dict(format="kayles-response-screen-v1", fingerprint=fingerprint,
                  commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  limits=dict(per_process_wall_s=2, total_wall_s=20, address_space_mib=512, threads=1, cpu=cpu),
                  status="running", runs=[], solver_processes=0)

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
        resource.setrlimit(resource.RLIMIT_CPU, (3,)*2)

    for mode, max_n in (("strict", args.reference_max_n), ("response", args.max_n)):
        remaining = 20-(time.monotonic()-started)
        assert remaining > 0, "batch budget reached"
        before = time.monotonic()
        command = [sys.executable, str(Path(__file__).resolve()), "--worker", mode,
                   "--max-n", str(max_n), "--pair-limit", str(args.pair_limit)]
        process = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   preexec_fn=limits, start_new_session=True)
        report["solver_processes"] += 1
        status = "complete"
        try:
            stdout, stderr = process.communicate(timeout=min(2, remaining))
            if process.returncode:
                status = "incomplete"
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
            status = "timeout"
        rows = [json.loads(line) for line in stdout.splitlines()]
        completed = [row for row in rows if row["status"] == "complete"]
        assert [row["n"] for row in completed] == list(range(1, len(completed)+1))
        for row in completed:
            if row["n"] in known:
                assert row["SG"] == known[row["n"]]
        assert all("SG" not in row for row in rows if row["status"] != "complete")
        if status == "complete":
            assert len(completed) == max_n
        run = dict(mode=mode, requested_max_n=max_n, status=status, returncode=process.returncode,
                   wall_s=round(time.monotonic()-before, 6), rows=rows, stderr=stderr)
        report["runs"].append(run)
        atomic_write_json(args.output, report)
        print(json.dumps(run), flush=True)
        if status != "complete":
            report["status"] = "candidate_stopped" if mode == "response" else "reference_stopped"
            break
    if report["status"] == "running":
        report["status"] = "complete"
    report["wall_s"] = round(time.monotonic()-started, 6)
    atomic_write_json(args.output, report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
