#!/usr/bin/env python3
"""Short proof-cost experiments; offline oracle or existing fast-search cache.

Each policy stops at its first failed size. The generator and independent
oracle-free checker both have 2 s wall/CPU/memory limits by default. Save only
complete checked certificates; no failed optimisation is reported as a proof.
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


def bounded(command, env, timeout):
    start = time.monotonic()
    child = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             env=env, start_new_session=True, preexec_fn=lambda: limits(512, timeout))
    try:
        stdout, stderr = child.communicate(timeout=timeout)
        status = "complete" if child.returncode == 0 else "failed"
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGKILL); stdout, stderr = child.communicate(); status = "timeout"
    return dict(status=status, returncode=child.returncode, wall_s=round(time.monotonic()-start, 6),
                stdout=stdout.strip(), stderr=stderr.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[8, 10, 12, 14, 16, 18])
    parser.add_argument("--engine", choices=("offline", "target", "lazy"), default="offline")
    parser.add_argument("--policies", nargs="+", choices=("tree_min", "first"), default=["tree_min", "first"])
    parser.add_argument("--per-run", type=float, default=2)
    parser.add_argument("--total", type=float, default=20)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--proof-dir", type=Path, required=True)
    args = parser.parse_args()
    known = known_values()
    if (args.sizes != sorted(set(args.sizes)) or not args.sizes or any(n not in known for n in args.sizes)
            or any(not math.isfinite(t) or t <= 0 for t in (args.per_run, args.total))):
        parser.error("increasing known widths and positive finite budgets required")
    cpu = min(os.sched_getaffinity(0)); os.sched_setaffinity(0, {cpu})
    source = {"offline": "src/cfast/strip3_proof_cost.c", "target": "src/cfast/strip3_target_proof.c",
              "lazy": "src/cfast/strip3_lazy_proof.c"}[args.engine]
    sources = list(dict.fromkeys([source, "src/cfast/strip3_proof_cost.c", "src/cfast/strip3_target.c",
               "src/cfast/strip3.c", "src/cfast/strip3_comb_tables.h",
               "tests/verify_strip3_equality_proof.py", "scripts/screen_strip3_proofs.py"]))
    if args.engine == "lazy": sources.extend(["src/cfast/strip3_lazy.c", "src/cfast/strip3_target_proof.c"])
    report = dict(format="strip3-offline-proof-screen-v1", status="running", samples=[],
                  configuration=dict(sizes=args.sizes, policies=args.policies, engine=args.engine),
                  limits=dict(per_process_s=args.per_run, total_s=args.total, address_space_mib=512,
                              cpu=cpu, threads=1, oracle_included=args.engine == "offline"),
                  interpretation="Full-SG-guided offline or available-proof tree optimisation; selected DAGs are not minimum DAGs.",
                  provenance=dict(commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                                  source_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources}))
    start = time.monotonic()
    def save():
        report["wall_s"] = round(time.monotonic()-start, 6)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n")
    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
    env.update(STRIP3_HT_BITS="22", STRIP3_HT_PREALLOC="1", STRIP3_THREADS="1", STRIP3_PROOF_BITS="20")
    if args.engine in ("target", "lazy"):
        env.update(STRIP3_TARGET_BITS="22", STRIP3_TARGET_FRONTIER="4", STRIP3_TARGET_CACHE_FIRST="1",
                   STRIP3_TARGET_SHARE_OPTIONS="1", STRIP3_TARGET_ORDER="2")
    if args.engine == "lazy": env["STRIP3_LAZY_QUANTUM"] = "256"
    report["configuration"]["environment"] = {k: v for k, v in env.items() if k.startswith("STRIP3_")}
    args.proof_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="strip3-proof-screen-") as folder:
        directory = Path(folder); binary = directory / "proof_cost"
        try:
            subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                            str(ROOT / sources[0])], check=True, capture_output=True, timeout=args.total)
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
            report.update(status="build_failed", error=str(error)); save(); return
        for policy in args.policies:
            for n in args.sizes:
                remaining = args.total-(time.monotonic()-start)
                if remaining <= 0: report["status"] = "total_budget"; save(); return
                cert = directory / "certificate.json"
                if cert.exists(): cert.unlink()
                row = bounded([str(binary), "board", str(n), "--certificate", str(cert)],
                              {**env, "STRIP3_PROOF_FIRST": str(int(policy == "first"))}, min(args.per_run, remaining))
                row.update(n=n, policy=policy)
                if row["status"] == "complete":
                    metrics = json.loads(row["stdout"].splitlines()[-1])
                    if metrics["n"] != n or metrics["SG"] != known[n]: row["status"] = "wrong_value"
                    else:
                        remaining = args.total-(time.monotonic()-start)
                        if remaining <= 0: row["status"] = "verification_budget"
                        else:
                            checked = bounded(["python3", str(ROOT / "tests/verify_strip3_equality_proof.py"), str(cert)],
                                              env, min(args.per_run, remaining))
                            row["verification"] = checked
                            if checked["status"] != "complete": row["status"] = "verification_"+checked["status"]
                            else:
                                summary = json.loads(checked["stdout"])
                                for key in ("n", "SG", "tree_cost", "sum_rules", "edges"):
                                    assert summary[key] == metrics[key], key
                                assert summary["queries"] == metrics["certificate_queries"]
                                raw = cert.read_bytes()
                                path = args.proof_dir / f"strip3_equality_{args.engine}_{policy}_{n}.json.gz"
                                path.write_bytes(gzip.compress(raw, mtime=0))
                                row.update(metrics=metrics, checked=summary, certificate=str(path),
                                           certificate_sha256=hashlib.sha256(raw).hexdigest(),
                                           certificate_bytes=len(raw), compressed_bytes=path.stat().st_size)
                report["samples"].append(row); save()
                print(json.dumps(dict(n=n, policy=policy, status=row["status"], wall_s=row["wall_s"],
                                      checked=row.get("checked"))), flush=True)
                if row["status"] != "complete": break
    report["status"] = "screened_requested_policies"; save()


if __name__ == "__main__":
    main()
