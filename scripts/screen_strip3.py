#!/usr/bin/env python3
"""Fast local screening of a board-N CLI against established exact SG values.

Hard wall/CPU/address-space limits; one CPU; fresh processes for uncached jobs;
stop at the first candidate failure. Checked completed jobs are reused unless
--fresh is supplied. Cache lookup time is not solver time or scaling evidence.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import signal
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
BUILD_FLAGS = ["-O3", "-march=native", "-Wall", "-Wextra", "-pthread"]
CACHE_FORMAT = "strip3-target100-cache-v1"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write_json(path, payload):
    """Write a report so an interrupted run cannot leave a partial cache."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(handle, "w") as stream:
            json.dump(payload, stream, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def relative_or_absolute(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def source_provenance(candidate_source, skip_reference):
    # Follow local includes recursively: several candidates include another .c
    # file which in turn includes strip3.c and its generated tables.
    pending = [candidate_source.resolve()]
    if not skip_reference:
        pending.append(ROOT / "src/cfast/strip3.c")
    paths = {Path(__file__).resolve()}
    while pending:
        path = pending.pop().resolve()
        if path in paths:
            continue
        paths.add(path)
        for include in re.findall(r'^\s*#\s*include\s*"([^"\n]+)"', path.read_text(), re.M):
            dependency = path.parent / include
            if not dependency.is_file():
                raise ValueError(f"cannot fingerprint local include {include} in {path}")
            pending.append(dependency)
    compiler = Path(shutil.which("gcc")).resolve()
    # CPU model/flags affect -march=native; MHz and other volatile fields do not.
    cpu = "\n".join(sorted(set(re.findall(r'^(?:model name|flags|Features)\s*:.*$',
                                            Path("/proc/cpuinfo").read_text(), re.M))))
    build_environment = {key: os.environ[key] for key in
                         ("CPATH", "C_INCLUDE_PATH", "LIBRARY_PATH", "COMPILER_PATH",
                          "GCC_EXEC_PREFIX", "LD_LIBRARY_PATH", "LD_PRELOAD", "SOURCE_DATE_EPOCH")
                         if key in os.environ}
    return {"candidate_source": relative_or_absolute(candidate_source.resolve()),
            "source_sha256": {relative_or_absolute(path): sha(path) for path in sorted(paths)},
            "build": {"compiler": str(compiler), "compiler_sha256": sha(compiler),
                      "flags": BUILD_FLAGS, "cpu_sha256": hashlib.sha256(cpu.encode()).hexdigest(),
                      "environment_sha256": hashlib.sha256(json.dumps(build_environment, sort_keys=True).encode()).hexdigest()},
            "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                               text=True).strip()}


def cache_fingerprint(provenance, configuration, memory_mib, per_run):
    """Append sizes without invalidating old rows; commits are provenance only."""
    return {"format": CACHE_FORMAT,
            "provenance": {key: value for key, value in provenance.items() if key != "commit"},
            "configuration": {key: value for key, value in configuration.items()
                               if key != "sizes"},
            "limits": {"address_space_mib": memory_mib, "threads": 1,
                       "per_process_wall_s": per_run}}


def load_cached_rows(path, fingerprint, known):
    """Return only checked complete rows compatible with this experiment."""
    if path is None:
        return {}, {"enabled": False}
    metadata = {"enabled": True, "path": str(path), "hits": 0}
    try:
        cached = json.loads(path.read_text())
    except FileNotFoundError:
        metadata["state"] = "missing"
        return {}, metadata
    except (OSError, ValueError) as error:
        metadata.update(state="invalid", reason=str(error))
        return {}, metadata
    if not isinstance(cached, dict) or cached.get("cache_fingerprint") != fingerprint:
        metadata.update(state="incompatible")
        return {}, metadata
    if not isinstance(cached.get("completed_samples"), list):
        metadata.update(state="invalid", reason="completed_samples must be an array")
        return {}, metadata
    rows = {}
    for row in cached["completed_samples"]:
        if not isinstance(row, dict) or row.get("status") != "complete":
            continue
        solver, n, sg = row.get("solver"), row.get("n"), row.get("SG")
        if solver not in ("candidate", "reference") or type(n) is not int or type(sg) is not int:
            continue
        if known.get(n) != sg or row.get("returncode") != 0:
            continue
        if str(row.get("stdout", "")).split() != ["3", str(n), str(sg)]:
            continue
        key = (solver, n)
        if key not in rows:
            checked = dict(row)
            rows[key] = checked
    metadata["state"] = "usable"
    metadata["available_rows"] = len(rows)
    return rows, metadata


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
    parser.add_argument("--cache", type=Path,
                        help="incremental JSON cache; defaults to --output")
    parser.add_argument("--fresh", action="store_true",
                        help="ignore existing completed rows but still write a new cache")
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
    cache_path = args.cache or args.output
    configuration = {"order": args.order, "low_first": args.low_first,
                     "sizes": args.sizes, "skip_reference": args.skip_reference,
                     "extra_env": extra_env}
    provenance = source_provenance(args.candidate_source, args.skip_reference)
    fingerprint = cache_fingerprint(provenance, configuration, args.memory_mib, args.per_run)
    if args.fresh:
        # A fresh request remeasures the requested keys, while preserving other
        # compatible completed sizes so a narrow refresh cannot erase progress.
        previous_rows, _ = load_cached_rows(cache_path, fingerprint, known)
        cached_rows, cache_info = {}, {"enabled": True, "path": str(cache_path),
                                      "state": "forced_fresh", "hits": 0}
        completed_rows = dict(previous_rows)
    else:
        cached_rows, cache_info = load_cached_rows(cache_path, fingerprint, known)
        completed_rows = dict(cached_rows)
    result = {"format": "strip3-target100-screen-v1", "status": "running", "samples": [],
              "configuration": configuration,
              "limits": {"per_process_wall_s": args.per_run, "total_wall_s": args.total,
                         "address_space_mib": args.memory_mib, "cpu": cpu, "threads": 1},
              "provenance": provenance,
              "cache_fingerprint": fingerprint,
              "cache": cache_info,
              "work": {"builds": 0, "solver_processes": 0},
              "interpretation": "One sample per size; finite measurements do not establish asymptotic complexity.",
              "timing_caveat": "Other local work may run concurrently; compare operation counts before timing ratios."}

    def save():
        result["wall_s"] = round(time.monotonic()-start, 6)
        # Preserve completed sizes outside this invocation's request, including
        # rows beyond its first failure. samples contains only current requests.
        result["completed_samples"] = sorted(completed_rows.values(),
                                             key=lambda row: (row["n"], row["solver"]))
        atomic_write_json(args.output, result)
        if cache_path.resolve() != args.output.resolve():
            atomic_write_json(cache_path, result)

    def record(row):
        result["samples"].append(row)
        if row["status"] == "complete":
            completed_rows[row["solver"], row["n"]] = row
    save()

    reference_active = not args.skip_reference
    needed = {name: any((name, n) not in cached_rows for n in args.sizes)
              for name in ("candidate", "reference")}
    with tempfile.TemporaryDirectory(prefix="strip3-screen-") as folder:
        binaries = {}
        sources = []
        if needed["candidate"]:
            sources.append(("candidate", args.candidate_source))
        if not args.skip_reference and needed["reference"]:
            sources.append(("reference", ROOT / "src/cfast/strip3.c"))
        for name, source in sources:
            binary = Path(folder) / name
            remaining = args.total - (time.monotonic()-start)
            if remaining <= 0:
                result["status"] = "total_budget"; save(); return
            try:
                result["work"]["builds"] += 1
                subprocess.run([provenance["build"]["compiler"], *BUILD_FLAGS, "-o", str(binary),
                                str(source)], check=True, capture_output=True, text=True, timeout=remaining)
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
                result["status"] = "build_failed"; result["error"] = str(error); save(); return
            binaries[name] = binary
        for n in args.sizes:
            for name in ("candidate", "reference"):
                if name == "reference" and not reference_active:
                    continue  # do not repeatedly spend the timeout on the old algorithm
                if (name, n) in cached_rows:
                    cached = dict(cached_rows[name, n], cached=True)
                    result["cache"]["hits"] += 1
                    record(cached)
                    print(json.dumps(cached), flush=True)
                    continue
                remaining = args.total - (time.monotonic()-start)
                if remaining <= 0:
                    result["status"] = "total_budget"; save(); return
                env = {"STRIP3_TARGET_BITS": "22", "STRIP3_TARGET_ORDER": str(args.order),
                       "STRIP3_TARGET_LOW_FIRST": str(args.low_first),
                       "STRIP3_THREADS": "1", "STRIP3_HT_BITS": "24", "STRIP3_HT_PREALLOC": "1",
                       "STRIP3_REDUCTIONS": "7", "STRIP3_PROGRESS": "0"}
                env.update(extra_env)
                timeout = min(args.per_run, remaining)
                result["work"]["solver_processes"] += 1
                row = measure(binaries[name], n, timeout, args.memory_mib, env)
                row["solver"] = name
                row["cached"] = False
                row["measured_at_commit"] = provenance["commit"]
                row["measurement_limits"] = dict(result["limits"], per_process_wall_s=timeout)
                if row["status"] == "complete" and row["SG"] != known[n]:
                    row["status"] = "wrong_value"
                record(row)
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
