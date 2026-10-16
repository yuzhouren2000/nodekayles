"""Fast-search certificates verify without exact-value calls during extraction."""
import json
import os
from pathlib import Path
import resource
import subprocess
import tempfile

from verify_strip3_equality_proof import verify

ROOT = Path(__file__).resolve().parents[1]


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (3,)*2)


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    known = {int(n): int(g) for _, n, g in map(str.split, (ROOT / "data/raw/strip3_3xn.tsv").read_text().splitlines())}
    with tempfile.TemporaryDirectory(prefix="strip3-target-proof-test-") as folder:
        directory = Path(folder); binary = directory / "target_proof"
        subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3_target_proof.c")], check=True, timeout=2)
        env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
        env.update(STRIP3_TARGET_BITS="18", STRIP3_PROOF_BITS="18")
        for mode in (0, 1):
            for n in range(1, 13):
                cert = directory / "proof.json"
                settings = dict(STRIP3_PROOF_FIRST=str(mode))
                if mode:
                    settings.update(STRIP3_TARGET_CACHE_FIRST="1", STRIP3_TARGET_SHARE_OPTIONS="1", STRIP3_TARGET_FRONTIER="4")
                result = subprocess.run([str(binary), "board", str(n), "--certificate", str(cert)], env={**env, **settings},
                                        text=True, capture_output=True, timeout=2, preexec_fn=limits)
                assert result.returncode == 0, result.stderr
                board, raw = result.stdout.splitlines()
                assert list(map(int, board.split())) == [3, n, known[n]]
                summary = json.loads(raw); assert summary["oracle_calls"] == 0
                checked = verify(json.loads(cert.read_text()))
                for key in ("n", "SG", "tree_cost", "sum_rules", "edges"): assert checked[key] == summary[key]
                assert checked["queries"] == summary["certificate_queries"]
        for extra in (dict(STRIP3_TARGET_BITS="8"), dict(STRIP3_PROOF_BITS="8")):
            cert = directory / "failed.json"
            result = subprocess.run([str(binary), "board", "12", "--certificate", str(cert)], env={**env, **extra},
                                    text=True, capture_output=True, timeout=2, preexec_fn=limits)
            assert result.returncode == 3 and not result.stdout and not cert.exists(), result
    print("24 exact-strip certificates verify independently; extraction performs zero exact-SG oracle calls; failure emits neither board nor proof")


if __name__ == "__main__":
    main()
