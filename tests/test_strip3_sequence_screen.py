"""Verify that a killed prefix preserves only completed, checked SG rows."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from screen_strip3_sequence import jobs, measure, references


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    known = references()
    with tempfile.TemporaryDirectory(prefix="strip3-sequence-test-") as folder:
        directory = Path(folder)
        labels, text = jobs("F", 2, known)
        for body, expected in (
                ("print(2, flush=True)\nwhile True: pass\n", "timeout"),
                ("print(999, flush=True)\n", "wrong_value"),
                ("print('bad', flush=True)\n", "invalid_output"),
                ("print(2, flush=True)\n", "incomplete_output")):
            fake = directory / "fake"
            fake.write_text("#!/usr/bin/env python3\n"+body); fake.chmod(0o700)
            row = measure(fake, labels, text, known, os.environ, .15)
            assert row["status"] == expected, row
            if expected == "timeout":
                assert row["completed"] == [dict(kind="F", n=1, SG=2, expected=2)]
                assert row["next_job"] == ["F", 2]
        binary = directory / "target"
        subprocess.run(["gcc", "-O3", "-pthread", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3_target.c")], check=True, timeout=2)
        labels, text = jobs("ADF", 4, known)
        row = measure(binary, labels, text, known,
                      {**os.environ, "STRIP3_TARGET_BATCH_PROGRESS": "1", "STRIP3_TARGET_BITS": "12"}, 2)
        assert row["status"] == "complete" and len(row["completed"]) == len(row["stages"]) == 12, row
        assert [s["memo"] for s in row["stages"]] == sorted(s["memo"] for s in row["stages"])
        invalid = subprocess.run([str(binary), "masks"], input="1 7\n", text=True, capture_output=True,
                                 env={**os.environ, "STRIP3_TARGET_BATCH_PROGRESS": "2"}, timeout=2)
        assert invalid.returncode == 2 and not invalid.stdout
    print("killed prefix keeps completed jobs only; wrong, invalid and incomplete outputs rejected; shared-memo progress matches all 12 jobs")


if __name__ == "__main__":
    main()
