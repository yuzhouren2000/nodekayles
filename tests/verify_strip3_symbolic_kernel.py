#!/usr/bin/env python3
"""Compile the standalone C kernel checker once and verify saved certificates.

No game solver is imported or run. Gzip is only a storage format. Each build
or check has a two-second deadline, one CPU and 512 MiB; the batch has 20 s.
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
    parser.add_argument("certificates", type=Path, nargs="+")
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    cpu = min(os.sched_getaffinity(0)); os.sched_setaffinity(0, {cpu})
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    started = time.monotonic()
    source = ROOT/"tests/verify_strip3_symbolic_kernel.c"
    report = {"format": "strip3-symbolic-kernel-audit-v1", "samples": [],
              "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (source, Path(__file__).resolve())},
              "limits": {"per_process_wall_s": 2, "total_wall_s": 20, "address_space_mib": 512,
                         "cpu": cpu, "threads": 1,
                         "total_budget_checks": "between stages; gzip reads and audit writes have no separate hard interrupt"}}
    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}

    def stage(command, text=None):
        remaining = 20-(time.monotonic()-started)
        if remaining <= 0:
            return {"status": "total_budget"}
        begin = time.monotonic()
        process = subprocess.Popen(list(map(str, command)), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True, env=env, start_new_session=True, preexec_fn=limits)
        try:
            stdout, stderr = process.communicate(text, timeout=min(2, remaining))
            status = "complete" if process.returncode == 0 else "failed"
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
            status = "timeout"
        return {"status": status, "returncode": process.returncode,
                "wall_s": round(time.monotonic()-begin, 6), "stdout": stdout.strip(), "stderr": stderr.strip()}

    with tempfile.TemporaryDirectory(prefix="strip3-kernel-check-") as folder:
        checker = Path(folder)/"check"
        report["build"] = stage(["gcc", "-O3", "-Wall", "-Wextra", "-Werror", "-o", checker, source])
        okay = report["build"]["status"] == "complete"
        if okay:
            for path in args.certificates:
                raw = path.read_bytes()
                data = gzip.decompress(raw) if path.suffix == ".gz" else raw
                result = stage([checker], data.decode("ascii"))
                row = {"path": str(path), "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data),
                       "file_sha256": hashlib.sha256(raw).hexdigest(), **result}
                if result["status"] == "complete":
                    row["checked"] = json.loads(result["stdout"])
                report["samples"].append(row)
                print(json.dumps(row), flush=True)
                if result["status"] != "complete":
                    okay = False; break
    report.update(status="verified" if okay else "incomplete", wall_s=round(time.monotonic()-started, 6))
    if args.audit:
        args.audit.parent.mkdir(parents=True, exist_ok=True)
        args.audit.write_text(json.dumps(report, indent=2)+"\n")
    if not okay:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
