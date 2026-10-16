#!/usr/bin/env python3
"""Short first-failure screen; values are checked only after a job completes.

Build and each run/check: one CPU, 512 MiB, two wall/CPU seconds. Stop the
whole batch on the first failure. Store completed certificates for reuse;
never repeat their search just to verify the recorded SG values.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (2,)*2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("availability", "history"), required=True)
    parser.add_argument("--generalise", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--proofs", type=Path, required=True)
    args = parser.parse_args()
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    started = time.monotonic()
    source = ROOT/"src/cfast"/f"strip3_{'history' if args.backend == 'history' else 'demand'}_symbolic.c"
    checker_source = ROOT/"tests/verify_strip3_symbolic_kernel.c"
    sources = (source, ROOT/"src/cfast/strip3_demand_symbolic.c", ROOT/"src/cfast/strip3_symbolic.c",
               checker_source, Path(__file__).resolve())
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=True).stdout.strip()
    report = {"format": "strip3-demand-screen-v1", "backend": args.backend, "generalise": args.generalise,
              "revision": revision, "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                                       for p in sources}, "samples": [],
              "limits": {"per_process_wall_s": 2, "per_process_cpu_s": 2, "total_wall_s": 20,
                         "address_space_mib": 512, "cpu": cpu, "threads": 1,
                         "total_budget_checks": "between stages; file I/O has no separate hard interrupt"}}
    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
    env.update(STRIP3_DEMAND_GENERALISE=str(args.generalise), STRIP3_SYMBOLIC_BITS="20",
               STRIP3_SYMBOLIC_GC="1", STRIP3_DEMAND_TRACE="1")

    def stage(command, data=None):
        remaining = 20-(time.monotonic()-started)
        if remaining <= 0:
            return {"status": "total_budget"}
        begin = time.monotonic()
        process = subprocess.Popen(list(map(str, command)), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True, env=env, start_new_session=True, preexec_fn=limits)
        try:
            stdout, stderr = process.communicate(data, timeout=min(2, remaining))
            status = "complete" if process.returncode == 0 else "failed"
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
            status = "timeout"
        return {"status": status, "returncode": process.returncode,
                "wall_s": round(time.monotonic()-begin, 6), "stdout": stdout.strip(), "stderr": stderr.strip()}

    okay = True
    with tempfile.TemporaryDirectory(prefix="strip3-demand-screen-") as folder:
        folder = Path(folder)
        binary, checker = folder/"solve", folder/"check"
        report["builds"] = []
        for output, code in ((binary, source), (checker, checker_source)):
            result = stage(["gcc", "-O3", "-Wall", "-Wextra", "-Werror", "-o", output, code])
            report["builds"].append(result)
            if result["status"] != "complete":
                okay = False
                break
        # One engine/configuration per batch. A failed job stops the batch.
        for n in ((6, 8, 10, 12, 14, 16) if okay else ()):
            proof = folder/"kernel.txt"
            row = {"n": n, "solve": stage([binary, "board", n, "--certificate", proof])}
            report["samples"].append(row)
            result = row["solve"]
            if result["status"] != "complete":
                okay = False
                print(json.dumps(row), flush=True)
                break
            parts = result["stdout"].split()
            assert len(parts) == 3 and parts[:2] == ["3", str(n)], result
            # The reference never enters the producer or changes its branches.
            references = {6: 3, 8: 2, 10: 3, 12: 5, 14: 4, 16: 3}
            assert int(parts[2]) == references[n], result
            data = proof.read_bytes()
            row["check"] = stage([checker], data.decode("ascii"))
            if row["check"]["status"] != "complete":
                okay = False
                print(json.dumps(row), flush=True)
                break
            row["checked"] = json.loads(row["check"]["stdout"])
            assert row["checked"]["SG"] == references[n] and row["checked"]["n"] == n
            args.proofs.mkdir(parents=True, exist_ok=True)
            saved = args.proofs/f"strip3_kernel_{n}.txt.gz"
            packed = gzip.compress(data, mtime=0)
            saved.write_bytes(packed)
            row["certificate"] = {"path": str(saved), "bytes": len(data), "gzip_bytes": len(packed),
                                  "sha256": hashlib.sha256(data).hexdigest(),
                                  "file_sha256": hashlib.sha256(packed).hexdigest()}
            print(json.dumps(row), flush=True)
    report.update(status="complete" if okay else "stopped_at_first_failure", wall_s=round(time.monotonic()-started, 6))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    main()
