#!/usr/bin/env python3
"""Bounded construction and standalone checking of symbolic SG kernels.

Every child gets one CPU, 512 MiB and two seconds by default. Stop at the
first incomplete construction or verification; never save an unchecked SG as
a completed result. The checker has its own BDD and domain construction.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time

from screen_strip3 import known_values, limits

ROOT = Path(__file__).resolve().parents[1]


def bounded(command, env, timeout, memory, stdin=None):
    start = time.monotonic()
    child = subprocess.Popen(list(map(str, command)), stdin=stdin, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True, env=env, start_new_session=True,
                             preexec_fn=lambda: limits(memory, timeout))
    try:
        stdout, stderr = child.communicate(timeout=timeout)
        status = "complete" if child.returncode == 0 else "failed"
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGKILL)
        stdout, stderr = child.communicate()
        status = "timeout"
    return {"status": status, "returncode": child.returncode,
            "wall_s": round(time.monotonic()-start, 6), "stdout": stdout.strip(), "stderr": stderr.strip()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[6, 8, 10, 12, 14])
    parser.add_argument("--engine", choices=("reachable", "partial"), default="reachable")
    parser.add_argument("--partial-initial", type=int, default=7)
    parser.add_argument("--bits", type=int, default=20)
    parser.add_argument("--checker-bits", type=int, default=20)
    parser.add_argument("--gc", type=int, choices=(0, 1), default=1)
    parser.add_argument("--per-run", type=float, default=2)
    parser.add_argument("--total", type=float, default=20)
    parser.add_argument("--memory-mib", type=int, default=512)
    parser.add_argument("--proof-directory", type=Path, default=ROOT/"data/proofs/strip3_symbolic_kernels")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    known = known_values()
    if (not args.sizes or args.sizes != sorted(set(args.sizes)) or any(n not in known for n in args.sizes)
            or not 8 <= args.bits <= 23 or not 8 <= args.checker_bits <= 22 or args.memory_mib < 128
            or not 0 <= args.partial_initial <= 300
            or not all(math.isfinite(v) and v > 0 for v in (args.per_run, args.total))):
        parser.error("require increasing known sizes, valid BDD budgets and positive finite resource limits")
    cpu = min(os.sched_getaffinity(0)); os.sched_setaffinity(0, {cpu})
    start = time.monotonic()
    sources = [ROOT/path for path in ("src/cfast/strip3_reachable_symbolic.c", "src/cfast/strip3_symbolic.c",
                                    "tests/verify_strip3_symbolic_kernel.c", "scripts/screen_strip3_symbolic_proofs.py",
                                    "scripts/screen_strip3.py", "src/cfast/strip3_partial_symbolic.c")]
    report = {"format": "strip3-symbolic-kernel-screen-v1", "status": "running", "samples": [],
              "provenance": {"commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                             "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}},
              "configuration": {"sizes": args.sizes, "bits": args.bits, "checker_bits": args.checker_bits,
                                "gc": args.gc, "engine": args.engine, "partial_initial": args.partial_initial},
              "limits": {"per_process_wall_s": args.per_run, "total_wall_s": args.total,
                         "address_space_mib": args.memory_mib, "cpu": cpu, "threads": 1,
                         "total_budget_checks": "between stages; compression and JSON writes have no separate hard interrupt"},
              "interpretation": "Exact certificates for fixed finite boards; no claim of a finite automaton valid for all lengths."}

    def save():
        report["wall_s"] = round(time.monotonic()-start, 6)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n")

    def stage(command, env, stdin=None):
        remaining = args.total-(time.monotonic()-start)
        if remaining <= 0:
            return {"status": "total_budget"}
        return bounded(command, env, min(args.per_run, remaining), args.memory_mib, stdin)

    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
    env.update(STRIP3_SYMBOLIC_BITS=str(args.bits), STRIP3_SYMBOLIC_GC=str(args.gc), STRIP3_KERNEL_BITS=str(args.checker_bits),
               STRIP3_PARTIAL_INITIAL=str(args.partial_initial))
    with tempfile.TemporaryDirectory(prefix="strip3-symbolic-kernels-") as folder:
        directory = Path(folder); producer, checker = directory/"producer", directory/"checker"
        for binary, source in ((producer, sources[0] if args.engine == "reachable" else sources[-1]), (checker, sources[2])):
            built = stage(["gcc", "-O3", "-Wall", "-Wextra", "-Werror", "-o", binary, source], env)
            if built["status"] != "complete":
                report["status"] = "build_failed"; report["build"] = built; save(); return
        for n in args.sizes:
            cert = directory/"kernel.txt"
            sample = {"n": n, "status": "incomplete"}
            sample["construction"] = result = stage([producer, "board", n, "--certificate", cert], env)
            if result["status"] == "complete":
                fields = result["stdout"].split()
                if fields != ["3", str(n), str(known[n])]:
                    sample["status"] = "invalid_claim"
                else:
                    with cert.open() as stream:
                        sample["verification"] = verified = stage([checker], env, stream)
                    if verified["status"] == "complete":
                        checked = json.loads(verified["stdout"])
                        if checked["n"] != n or checked["SG"] != known[n]:
                            sample["status"] = "checker_disagreement"
                        else:
                            data = cert.read_bytes()
                            args.proof_directory.mkdir(parents=True, exist_ok=True)
                            prefix = "strip3_kernel" if args.engine == "reachable" else "strip3_partial_kernel"
                            destination = args.proof_directory/f"{prefix}_{n}.txt.gz"
                            compressed = gzip.compress(data, mtime=0)
                            destination.write_bytes(compressed)
                            sample.update(status="verified", SG=checked["SG"], checked=checked,
                                          certificate={"path": str(destination.resolve().relative_to(ROOT)),
                                                       "bytes": len(data), "gzip_bytes": len(compressed),
                                                       "sha256": hashlib.sha256(data).hexdigest(),
                                                       "gzip_sha256": hashlib.sha256(compressed).hexdigest()})
                    else:
                        sample["status"] = "verification_incomplete"
            report["samples"].append(sample); save()
            print(json.dumps(sample), flush=True)
            if sample["status"] != "verified":
                report["status"] = "screen_stopped"; save(); return
    report["status"] = "completed_requested_sizes"; save()


if __name__ == "__main__":
    main()
