"""Check exact history-image language, symbolic SG and BDD compaction.

The reference enumerates independent vertex sets and uses direct subset mex;
it shares neither the column transducer nor any BDD operation with the solver.
"""
import os
import json
from pathlib import Path
import resource
import subprocess
import tempfile

from test_strip3 import neighbourhoods, subset_table

ROOT = Path(__file__).resolve().parents[1]


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CPU, (2,)*2)


def run(binary, mode, n, *arguments, **settings):
    env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
    env.update(settings)
    return subprocess.run([str(binary), mode, str(n), *map(str, arguments)], capture_output=True, text=True,
                          timeout=2, env=env, preexec_fn=limits)


def independent_history_images(n):
    closed = neighbourhoods(n)
    full = (1 << (3*n))-1
    images = set()

    def visit(candidates, available):
        images.add(available)
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            v = bit.bit_length()-1
            visit(candidates & ~closed[v], available & ~closed[v])

    visit(full, full)
    return images


def sparse_oracle_certificate(n, images, table):
    """Small diagnostic only: use the independent SG table to select replies.

    This is not a new search algorithm or a speed comparison. It demonstrates
    that the certificate conditions accept partial lower kernels.
    """
    closed = neighbourhoods(n)
    full = (1 << (3*n))-1
    target = table[full]
    facts = [set() for _ in range(target+1)]

    def losing(mask, k):
        assert table[mask] == k
        if mask in facts[k]:
            return
        facts[k].add(mask)
        for v in range(3*n):
            if mask >> v & 1:
                winning(mask & ~closed[v], k)
        for j in range(k):
            winning(mask, j)

    def winning(mask, k):
        assert table[mask] != k
        if table[mask] < k:
            losing(mask, table[mask])
        else:
            child = next(mask & ~closed[v] for v in range(3*n)
                         if mask >> v & 1 and table[mask & ~closed[v]] == k)
            losing(child, k)

    losing(full, target)
    nodes, intern = [], {}

    def make(v, low, high):
        if low == high:
            return low
        flip = low & 1
        key = v, low ^ flip, high ^ flip
        if key not in intern:
            nodes.append(key)
            intern[key] = 2*len(nodes)
        return intern[key] ^ flip

    def build(members, v=0):
        if not members:
            return 0
        if v == 3*n:
            return 1
        return make(v, build([m for m in members if not (m >> v & 1)], v+1),
                    build([m for m in members if m >> v & 1], v+1))

    domain = build(sorted(images))
    roots = [build(sorted(group)) for group in facts]
    text = (f"node-kayles-symbolic-kernel-v1\n{n} {target} {len(nodes)+1}\n{domain}\n"
            + " ".join(map(str, roots))+"\n"
            + "".join(f"{v} {low} {high}\n" for v, low, high in nodes))
    assert any(len(facts[j]) < sum(table[m] == j for m in images) for j in range(target))
    return text, [len(group) for group in facts]


def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    expected_counts = (4, 14, 54, 172, 600, 2037)
    checked_values, checked_moves = 0, 0
    with tempfile.TemporaryDirectory(prefix="strip3-reachable-symbolic-test-") as folder:
        binary, baseline, checker, partial = [Path(folder)/name for name in ("reachable", "induced", "check", "partial")]
        for path, source in ((binary, ROOT/"src/cfast/strip3_reachable_symbolic.c"),
                             (baseline, ROOT/"src/cfast/strip3_symbolic.c"),
                             (checker, ROOT/"tests/verify_strip3_symbolic_kernel.c"),
                             (partial, ROOT/"src/cfast/strip3_partial_symbolic.c")):
            subprocess.run(["gcc", "-O3", "-Wall", "-Wextra", "-Werror", "-o", str(path),
                            str(source)], check=True, timeout=2, preexec_fn=limits)

        def check(text, **settings):
            env = {k: v for k, v in os.environ.items() if not k.startswith("STRIP3_")}
            return subprocess.run([str(checker)], input=text, capture_output=True, text=True,
                                  timeout=2, preexec_fn=limits, env={**env, **settings})

        certificate_checks = 0
        for n in range(1, 7):
            images = independent_history_images(n)
            assert len(images) == expected_counts[n-1]
            result = run(binary, "domain", n)
            assert result.returncode == 0, result.stderr
            assert list(map(int, result.stdout.split())) == [int(m in images) for m in range(1 << (3*n))]
            # Every real successor remains in the independent-history image.
            closed = neighbourhoods(n)
            for mask in images:
                for v in range(3*n):
                    if mask >> v & 1:
                        assert mask & ~closed[v] in images
                        checked_moves += 1
            expected = subset_table(n)
            for gc, bits in (("0", "20"), ("1", "20"), ("1", "18")):
                result = run(binary, "table", n, STRIP3_SYMBOLIC_GC=gc, STRIP3_SYMBOLIC_BITS=bits)
                assert result.returncode == 0, result.stderr
                assert list(map(int, result.stdout.split())) == [g if m in images else -1
                                                                 for m, g in enumerate(expected)]
                checked_values += len(images)
                if n == 6 and gc == "1" and bits == "18":
                    assert "gc=0 " not in result.stderr, result.stderr
            for initial in ("0", "7"):
                result = run(partial, "table", n, STRIP3_PARTIAL_INITIAL=initial, STRIP3_SYMBOLIC_BITS="18")
                assert result.returncode == 0, result.stderr
                assert list(map(int, result.stdout.split())) == [g if m in images else -1
                                                                 for m, g in enumerate(expected)]
                checked_values += len(images)
            if n == 5:
                result = run(baseline, "table", n)
                assert result.returncode == 0, result.stderr
                assert list(map(int, result.stdout.split())) == expected
            cert = Path(folder)/"kernel.txt"
            result = run(binary, "board", n, "--certificate", cert)
            assert result.returncode == 0, result.stderr
            certificate_text = cert.read_text()
            verified = check(certificate_text)
            assert verified.returncode == 0, verified.stderr
            metrics = json.loads(verified.stdout)
            assert metrics["n"] == n and metrics["SG"] == expected[-1]
            assert metrics["predecessor_calls"] in (expected[-1]+1, expected[-1]+2)
            certificate_checks += 1
        for gc, bits in (("0", "20"), ("1", "18")):
            result = run(binary, "board", 8, "--certificate", cert,
                         STRIP3_SYMBOLIC_GC=gc, STRIP3_SYMBOLIC_BITS=bits)
            assert result.returncode == 0, result.stderr
            verified = check(cert.read_text())
            assert verified.returncode == 0 and json.loads(verified.stdout)["SG"] == 2, verified.stderr
            certificate_checks += 1
        for initial in ("0", "7"):
            result = run(partial, "board", 8, "--certificate", cert,
                         STRIP3_PARTIAL_INITIAL=initial, STRIP3_SYMBOLIC_BITS="18")
            assert result.returncode == 0, result.stderr
            verified = check(cert.read_text())
            assert verified.returncode == 0 and json.loads(verified.stdout)["SG"] == 2, verified.stderr
            certificate_checks += 1
        sparse_text, sparse_sizes = sparse_oracle_certificate(6, images, expected)
        verified = check(sparse_text)
        assert verified.returncode == 0, verified.stderr
        assert json.loads(verified.stdout)["SG"] == 3
        assert json.loads(verified.stdout)["complete_lower_kernels"] < 3
        certificate_checks += 1
        # Preserve the n=6 proof for semantic mutations with several lower kernels.
        result = run(binary, "board", 6, "--certificate", cert)
        assert result.returncode == 0, result.stderr
        original = cert.read_text().splitlines()
        mutations = []
        changed = original[:]; changed[2] = "1"; mutations.append(changed)  # wrong domain
        changed = original[:]; roots = changed[3].split(); roots[0] = "0"
        changed[3] = " ".join(roots); mutations.append(changed)  # incomplete lower kernel
        changed = original[:]; roots = changed[3].split(); roots[-1] = "0"
        changed[3] = " ".join(roots); mutations.append(changed)  # full board omitted
        changed = original[:]; roots = changed[3].split(); roots[-1] = changed[2]
        changed[3] = " ".join(roots); mutations.append(changed)  # invalid kernel
        # A singleton full-board set has no move to itself. The last-kernel
        # witness condition must still reject its false claim of SG=0.
        changed = original[:]
        n, _, count = map(int, changed[1].split())
        child = 1
        for v in reversed(range(3*n)):
            changed.append(f"{v} 0 {child}")
            child = 2*count; count += 1
        changed[1] = f"{n} 0 {count}"; changed[3] = str(child)
        bad = check("\n".join(changed)+"\n")
        assert bad.returncode != 0 and "outside its winning witnesses" in bad.stderr, bad
        for changed in mutations:
            bad = check("\n".join(changed)+"\n")
            assert bad.returncode != 0 and not bad.stdout
        for text in ("\n".join(original[:-1])+"\n", cert.read_text()+"trailing\n"):
            assert check(text).returncode != 0
        assert check(cert.read_text(), STRIP3_KERNEL_BITS="8").returncode != 0
        result = run(binary, "board", 4, "--certificate", Path(folder)/"absent"/"proof.txt")
        assert result.returncode == 3 and not result.stdout and "write failed" in result.stderr
        for n in (0, 101, "6bad"):
            assert run(binary, "board", n).returncode == 2
        for mode in ("domain", "table"):
            assert run(binary, mode, 7).returncode == 2
        for setting in (dict(STRIP3_SYMBOLIC_BITS="7"), dict(STRIP3_SYMBOLIC_GC="2")):
            result = run(binary, "board", 4, **setting)
            assert result.returncode == 2 and not result.stdout
        for gc in ("0", "1"):
            result = run(binary, "board", 12, STRIP3_SYMBOLIC_BITS="8", STRIP3_SYMBOLIC_GC=gc)
            assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr
            result = run(partial, "board", 12, STRIP3_SYMBOLIC_BITS="8", STRIP3_SYMBOLIC_GC=gc)
            assert result.returncode == 3 and not result.stdout and "incomplete" in result.stderr
        for initial in ("-1", "301", "4bad"):
            result = run(partial, "board", 4, STRIP3_PARTIAL_INITIAL=initial)
            assert result.returncode == 2 and not result.stdout
        assert run(baseline, "domain", 2).returncode == 2
    print({"domain_counts": expected_counts, "independent_SG_comparisons": checked_values,
           "real_move_closure_checks": checked_moves, "baseline_induced_SG_comparisons": 32768,
           "standalone_kernel_certificates": certificate_checks,
           "sparse_n6_kernel_facts": sparse_sizes,
           "semantic_tampering_and_missing_successor_witness": "rejected",
           "forced_compaction_and_budget_failures": "passed"})


if __name__ == "__main__":
    main()
