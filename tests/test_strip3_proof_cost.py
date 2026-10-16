"""Direct small-game minimum-tree recurrence and certificate mutation tests."""
import copy
from functools import cache
import json
import os
from pathlib import Path
import resource
import subprocess
import tempfile

from test_strip3 import subset_table
from verify_strip3_equality_proof import (verify, normal_parts, options, cells,
                                         axiom_value)

ROOT = Path(__file__).resolve().parents[1]


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (3,)*2)


def run(binary, args, text=None, **extra):
    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
    env.update(STRIP3_HT_BITS="20", STRIP3_PROOF_BITS="18")
    env.update(extra)
    return subprocess.run([str(binary), *args], input=text, text=True, capture_output=True,
                          timeout=2, preexec_fn=limits, env=env)


def minimum_reference(values):
    @cache
    def primitive(mask):
        for rule in (4, 2, 3):
            try: return axiom_value(mask, rule)
            except AssertionError: pass
        return None

    @cache
    def sum_cost(parts, k):
        total, offset = 1, 0
        for part in parts[1:]:
            offset ^= values[part]
            total += equality_cost(part, values[part])
        if parts: total += equality_cost(parts[0], k^offset)
        return total

    @cache
    def equality_cost(mask, k):
        if k > mask.bit_count() or primitive(mask) is not None: return 1
        g = values[mask]
        successor_sums = options(mask)
        if g == k:
            return 1+sum(sum_cost(parts, k) for parts in successor_sums)+sum(equality_cost(mask, j) for j in range(k))
        choices = []
        for parts in successor_sums:
            value = 0
            for part in parts: value ^= values[part]
            if value == k: choices.append(sum_cost(parts, k))
        if g < k: choices.append(equality_cost(mask, g))
        return 1+min(choices)

    return sum_cost


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    values = subset_table(3)
    reference = minimum_reference(values)
    jobs = [(mask, k) for mask in range(512) for k in range(7)]
    text = "".join(f"3 {mask:x} {k}\n" for mask, k in jobs)
    with tempfile.TemporaryDirectory(prefix="strip3-proof-test-") as folder:
        directory = Path(folder); binary = directory / "proof_cost"
        subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-pthread", "-o", str(binary),
                        str(ROOT / "src/cfast/strip3_proof_cost.c")], check=True, timeout=2)
        for first in (0, 1):
            result = run(binary, ["queries"], text, STRIP3_PROOF_FIRST=str(first))
            assert result.returncode == 0, result.stderr
            rows = list(map(json.loads, result.stdout.splitlines())); assert len(rows) == len(jobs)
            for (mask, k), row in zip(jobs, rows):
                assert row["g"] == values[mask] and row["equal"] is (values[mask] == k)
                expected = reference(normal_parts(mask), k)
                assert int(row["cost"]) >= expected
                if not first: assert int(row["cost"]) == expected, (mask, k, row, expected)
        for first in (0, 1):
            cert = directory / f"cert_{first}.json"
            result = run(binary, ["board", "8", "--certificate", str(cert)], STRIP3_PROOF_FIRST=str(first))
            assert result.returncode == 0, result.stderr
            report = json.loads(cert.read_text()); summary = json.loads(result.stdout)
            checked = verify(report)
            for key in ("n", "SG", "tree_cost", "sum_rules", "edges"): assert checked[key] == summary[key]
            assert checked["queries"] == summary["certificate_queries"]
            if first: continue
            # Reject omitted moves, altered offsets, invalid witnesses and a
            # false root claim. Each change attacks a proof obligation.
            all_node = next(r for r in report["nodes"] if r["rule"] == 5 and r["options"])
            mutations = []
            bad = copy.deepcopy(report); bad["nodes"][all_node["id"]-1]["options"].pop(); mutations.append(bad)
            bad = copy.deepcopy(report); bad["nodes"][0]["k"] += 1; mutations.append(bad)
            move = next(r for r in report["nodes"] if r["rule"] == 6)
            bad = copy.deepcopy(report); bad["nodes"][move["id"]-1]["move"] = 125; mutations.append(bad)
            record = next(r for r in report["nodes"] if r["rule"] == 6 and len(r["witness"]["parts"]) > 1)
            bad = copy.deepcopy(report); bad["nodes"][record["id"]-1]["witness"]["k"] ^= 1; mutations.append(bad)
            for bad in mutations:
                try: verify(bad)
                except (AssertionError, KeyError): pass
                else: raise AssertionError("invalid certificate accepted")
        for extra in (dict(STRIP3_PROOF_BITS="8"), dict(STRIP3_HT_BITS="8")):
            cert = directory / "failed.json"
            result = run(binary, ["board", "12", "--certificate", str(cert)], **extra)
            assert result.returncode == 3 and not result.stdout and not cert.exists(), result
        assert run(binary, ["queries"], "3 200 0\n").returncode == 2
    print("3584 small sum queries: all truth values and every minimum tree cost match independent subset/coordinate recurrence")
    print("both policies export oracle-free verified SG proofs; four proof corruptions and exhausted budgets are rejected")


if __name__ == "__main__":
    main()
