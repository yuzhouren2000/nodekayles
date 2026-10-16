"""Check that bounded strip screens are safely resumable."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from screen_strip3 import atomic_write_json, cache_fingerprint, load_cached_rows


def check_cache_validation():
    known = {12: 5, 16: 3}
    provenance = {"commit": "abc", "candidate_source": "fake.c",
                  "source_sha256": {"fake.c": "source-v1"}}
    configuration = {"order": 2, "low_first": 0, "sizes": [12, 16],
                    "skip_reference": True, "extra_env": {}}
    fingerprint = cache_fingerprint(provenance, configuration, 512, 2)
    assert fingerprint == cache_fingerprint(dict(provenance, commit="def"),
                                            dict(configuration, sizes=[12, 16, 18]), 512, 2)
    with tempfile.TemporaryDirectory(prefix="strip3-cache-test-") as folder:
        path = Path(folder) / "screen.json"
        atomic_write_json(path, {"cache_fingerprint": fingerprint, "completed_samples": [
            {"solver": "candidate", "n": 12, "status": "complete", "SG": 5,
             "returncode": 0, "stdout": "3 12 5"},
            {"solver": "candidate", "n": 16, "status": "timeout", "SG": 3,
             "returncode": -9, "stdout": ""},
            {"solver": "candidate", "n": 16, "status": "complete", "SG": 99,
             "returncode": 0, "stdout": "3 16 99"},
            {"solver": "candidate", "n": 16, "status": "complete", "SG": 3,
             "returncode": 0, "stdout": "3 16 3"},
            {"solver": "reference", "n": 12, "status": "complete", "SG": 5,
             "returncode": 1, "stdout": "3 12 5"},
            {"solver": "reference", "n": 12, "status": "complete", "SG": 5,
             "returncode": 0, "stdout": "3 12 5 extra"},
            {"solver": "reference", "n": 12, "status": "complete", "SG": 99,
             "returncode": 0, "stdout": "3 12 99"},
            {"solver": "reference", "n": 16, "status": "timeout", "SG": 3,
             "returncode": 0, "stdout": "3 16 3"},
            {"solver": "reference", "n": "12", "status": "complete", "SG": 5,
             "returncode": 0, "stdout": "3 12 5"},
            None,
        ]})
        rows, info = load_cached_rows(path, fingerprint, known)
        assert info["state"] == "usable" and info["available_rows"] == 2, info
        assert set(rows) == {("candidate", 12), ("candidate", 16)}, rows
        assert rows["candidate", 12]["SG"] == 5
        assert rows["candidate", 12]["status"] == "complete"
        assert rows["candidate", 16]["SG"] == 3
        assert rows["candidate", 16]["stdout"] == "3 16 3"

        for changed in (
                cache_fingerprint(dict(provenance, source_sha256={"fake.c": "v2"}),
                                  configuration, 512, 2),
                cache_fingerprint(provenance, dict(configuration, order=1), 512, 2),
                cache_fingerprint(provenance, dict(configuration, extra_env={"STRIP3_X": "1"}), 512, 2),
                cache_fingerprint(provenance, configuration, 256, 2),
                cache_fingerprint(provenance, configuration, 512, 1)):
            rows, info = load_cached_rows(path, changed, known)
            assert not rows and info["state"] == "incompatible", info

        rows, _ = load_cached_rows(path, fingerprint, {12: 0, 16: 3})
        assert set(rows) == {("candidate", 16)}, rows

        broken = Path(folder) / "broken.json"
        broken.write_text("not json")
        rows, info = load_cached_rows(broken, fingerprint, known)
        assert not rows and info["state"] == "invalid", info


def check_cli_resume():
    # Exercise real compiler/solver invocations, including a request that shrinks
    # and then grows again. The small positions keep this a seconds-only check.
    with tempfile.TemporaryDirectory(prefix="strip3-cache-cli-") as folder:
        output, cache = Path(folder) / "report.json", Path(folder) / "cache.json"

        def run(sizes, *extra):
            subprocess.run([sys.executable, str(ROOT / "scripts/screen_strip3.py"),
                            "--skip-reference", "--sizes", *map(str, sizes),
                            "--output", str(output), "--cache", str(cache), *extra],
                           check=True, capture_output=True, text=True, timeout=30)
            report = json.loads(output.read_text())
            assert report["status"] == "completed_requested_sizes", report
            assert report == json.loads(cache.read_text())
            assert [row["n"] for row in report["samples"]] == sizes
            return report

        report = run([12, 16])
        assert report["work"] == {"builds": 1, "solver_processes": 2}, report
        assert not any(row["cached"] for row in report["samples"])
        report = run([12, 16])
        assert report["work"] == {"builds": 0, "solver_processes": 0}, report
        assert report["cache"]["hits"] == 2
        assert all(row["cached"] for row in report["samples"])

        report = run([12, 16, 18])
        assert report["work"] == {"builds": 1, "solver_processes": 1}, report
        assert [row["cached"] for row in report["samples"]] == [True, True, False]
        report = run([12])
        assert report["work"] == {"builds": 0, "solver_processes": 0}, report
        assert len(report["completed_samples"]) == 3
        report = run([12, 16, 18])
        assert report["work"] == {"builds": 0, "solver_processes": 0}, report
        assert report["cache"]["hits"] == 3

        report = run([12], "--fresh")
        assert report["cache"]["state"] == "forced_fresh"
        assert report["work"] == {"builds": 1, "solver_processes": 1}, report
        assert len(report["completed_samples"]) == 3
        report = run([12], "--order", "1")
        assert report["cache"]["state"] == "incompatible"
        assert report["work"] == {"builds": 1, "solver_processes": 1}, report


def main():
    check_cache_validation()
    check_cli_resume()
    print("screen cache rejects invalid/stale rows; real CLI resumes, extends, shrinks and refreshes correctly")


if __name__ == "__main__":
    main()
