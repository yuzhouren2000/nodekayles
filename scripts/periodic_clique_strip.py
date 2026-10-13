#!/usr/bin/env python3
"""Exact defect-parameter interval algebra for any three-row induced strip.

Clique-column tails may visit all three rows; every move empties its column.
Their simultaneous finite period window lifts to infinity, when it exists.
Nonseparating moves in general strips are tracked per special column. See
notes/16_clique_column_strips.md for the proofs and precise scope.
"""

import argparse
from functools import cache
import hashlib
import json
from pathlib import Path

from strip3_defect_algebra import RootedComb, mex


def validate_word(word):
    word = tuple(word)
    if not word or any(type(c) is not int or c not in (0, 1, 2, 3, 4, 6) for c in word):
        raise ValueError("each column must be empty, a vertex, or a vertical edge")
    return word


class CliqueColumnIntervals:
    """Any immutable (not necessarily periodic) clique-column word."""

    def __init__(self, columns):
        self.columns = tuple(columns)
        if self.columns:
            validate_word(self.columns)
        self.value = cache(self._value)

    def _value(self, left, right, delete_left=0, delete_right=0):
        if left >= right:
            return 0
        options = set()
        for i in range(left, right):
            playable = self.columns[i]
            if i == left:
                playable &= ~delete_left
            if i == right - 1:
                playable &= ~delete_right
            for bit in (1, 2, 4):
                if playable & bit:
                    options.add(self.value(left, i, delete_left, bit) ^
                                self.value(i + 1, right, bit, delete_right))
        return mex(options)


class DefectiveThreeRowIntervals:
    """Exact FPT interval algebra for any three-row induced column word.

    Only columns 5 (nonadjacent outer vertices) and 7 (all three vertices)
    can support a move that does not empty its column. At each such original
    column this can happen at most once. The choice is unplayed/top/bottom.
    Endpoint deletions summarize all separating moves within an interval.
    """

    def __init__(self, columns):
        self.columns = tuple(columns)
        if any(type(c) is not int or not 0 <= c <= 7 for c in self.columns):
            raise ValueError("three-row column masks must be in 0..7")
        self.special = {i: j for j, i in enumerate(i for i, c in enumerate(self.columns)
                                                 if c in (5, 7))}
        self.value = cache(self._value)
        self.transformed = cache(self._transformed)

    def _transformed(self, choices):
        columns = list(self.columns)
        for i, j in self.special.items():
            choice = choices // (3 ** j) % 3
            if not choice:
                continue
            bit = 1 if choice == 1 else 4
            columns[i] &= ~(bit | (bit << 1) | (bit >> 1))
            if i:
                columns[i - 1] &= ~bit
            if i + 1 < len(columns):
                columns[i + 1] &= ~bit
        return tuple(columns)

    def _value(self, choices, left, right, delete_left=0, delete_right=0):
        if left >= right:
            return 0
        columns = self.transformed(choices)
        options = set()
        for i in range(left, right):
            playable = columns[i]
            if i == left:
                playable &= ~delete_left
            if i == right - 1:
                playable &= ~delete_right
            for bit in (1, 2, 4):
                if not playable & bit:
                    continue
                if not playable & ~(bit | (bit << 1) | (bit >> 1)):
                    options.add(self.value(choices, left, i, delete_left, bit) ^
                                self.value(choices, i + 1, right, bit, delete_right))
                else:
                    j = self.special[i]
                    assert choices // (3 ** j) % 3 == 0 and bit in (1, 4)
                    added = (1 if bit == 1 else 2) * 3 ** j
                    options.add(self.value(choices + added, left, right,
                                           delete_left, delete_right))
        return mex(options)

    def full_value(self):
        return self.value(0, 0, len(self.columns))


def ladder_endpoint_value(n, left=0, right=0):
    """Exact full two-row ladder table, including its unbounded corner defect."""
    if type(n) is not int or n < 0 or left not in range(4) or right not in range(4):
        raise ValueError("invalid ladder interval")
    if n == 0:
        return 0
    if n == 1:
        return int(bool(3 & ~(left | right)))
    if left == 3:
        return ladder_endpoint_value(n - 1, 0, right)
    if right == 3:
        return ladder_endpoint_value(n - 1, left, 0)
    if left == right == 0:
        return n % 2
    if left == 0 or right == 0:
        return n
    return int(left == right)


def window_end(start, period):
    return 2 * start + 2 * period - 1


def period_start(rows, period, minimum=1):
    """rows[n] contains EVERY state of that length, in fixed order."""
    available = len(rows) - 1
    if period < 1 or minimum < 1:
        raise ValueError("positive period and onset required")
    last_mismatch = minimum - 1
    for n in range(minimum, available - period + 1):
        if rows[n] != rows[n + period]:
            last_mismatch = n
    start = last_mismatch + 1
    return start if window_end(start, period) <= available else None


class PeriodicStrip:
    """All phases and arbitrary endpoint deletion masks are retained."""

    def __init__(self, word):
        self.word = validate_word(word)
        self.degree = len(self.word)
        self.rows = [[0] * (self.degree * 64)]

    def value(self, n, phase=0, left=0, right=0):
        return self.rows[n][(phase % self.degree) * 64 + left * 8 + right]

    def extend(self, maximum):
        for n in range(len(self.rows), maximum + 1):
            row = []
            for phase in range(self.degree):
                # Identical effective endpoint deletions share their result.
                endpoints = {}
                for left in range(8):
                    for right in range(8):
                        key = (left & self.word[phase],
                               right & self.word[(phase + n - 1) % self.degree])
                        if key not in endpoints:
                            endpoints[key] = mex({value for _, value in
                                                  self.options(n, phase, left, right)})
                        row.append(endpoints[key])
            self.rows.append(row)
        return self

    def options(self, n, phase=0, left=0, right=0, certificate=None):
        if certificate:
            start, period = certificate
            lookup = lambda length, p, l, r: self.lookup(length, p, l, r, start, period)
            positions = representative_columns(n, start, period)
        else:
            lookup = self.value
            positions = range(n)
        result = []
        for i in positions:
            playable = self.word[(phase + i) % self.degree]
            if i == 0:
                playable &= ~left
            if i == n - 1:
                playable &= ~right
            for bit in (1, 2, 4):
                if playable & bit:
                    value = lookup(i, phase, left, bit) ^ lookup(
                        n - i - 1, phase + i + 1, bit, right)
                    result.append(((bit.bit_length(), i + 1), value))
        return result

    def certificate(self, period, minimum=1):
        if period % self.degree:
            raise ValueError("period must preserve the column-word phase")
        return period_start(self.rows, period, minimum)

    def lookup(self, n, phase, left, right, start, period):
        index = n if n < start else start + (n - start) % period
        return self.value(index, phase, left, right)

    def report(self, period, start):
        end = window_end(start, period)
        rows = self.rows[:end + 1]
        if (period % self.degree or end >= len(self.rows)
                or any(self.rows[n] != self.rows[n + period]
                       for n in range(start, 2 * start + period))):
            raise ValueError("period window does not pass")
        return {"word": list(self.word), "period": period, "start": start,
                "maximum_length": end, "rows": rows,
                "rows_sha256": table_hash(rows),
                "full_period_blocks": [[self.value(n, phase) for n in
                                        range(start, start + period)]
                                       for phase in range(self.degree)]}


def table_hash(rows):
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()


def representative_columns(n, start, period):
    # Retain the short right factors as actual endpoint moves. All remaining
    # moves can shift left by period until their left factor is <start+period.
    return sorted(set(range(min(n, start + period))) |
                  set(range(max(0, n - start), n)))


def example_report():
    from strip3_defect_algebra import full_head
    strips, cores = [], []
    for word in ("23", "223", "23233", "4646313136"):
        period = 2 * len(word)
        if word == "4646313136":
            period = 10
        strip = PeriodicStrip(map(int, word)).extend(80)
        start = strip.certificate(period)
        if start is None:
            raise AssertionError("example lost its simultaneous period window")
        # The pure-tail certificate need not repeat a larger core window.
        strips.append(strip.report(period, start))
        if word == "4646313136":
            for head in (1, 2, 3, 4):
                vertices, adjacency = full_head(head)
                core = RootedPeriodicStrip(adjacency, strip).extend(80)
                onset = core.certificate(period)
                if onset is None:
                    raise AssertionError("rooted example lost its period window")
                cores.append({**core.report(period, onset), "word": list(strip.word),
                              "head": head, "cells": vertices,
                              "reachable_length_residues": [1, 2, 3, 4, 6, 7, 8, 9]})
    return {"format": "nodekayles-clique-strip-period-suite-v1",
            "strips": strips, "cores": cores}


class RootedPeriodicStrip:
    """Finite core joined by ONE edge to a vertex of the first tail column.

    Prefixes retain every right-end deletion and every root-preserving core
    descendant. No claim of unconditional stabilization is made here.
    """

    def __init__(self, adjacency, strip, attachment=None):
        if attachment is None:
            attachment = strip.word[0] & -strip.word[0]
        if attachment not in (1, 2, 4) or not strip.word[0] & attachment:
            raise ValueError("attachment vertex must exist in the first tail column")
        self.strip = strip
        self.core = RootedComb(adjacency)
        self.first = attachment
        self.rows = {mask: [[self.core.core_value(mask)] * 8] for mask in self.core.states}

    def value(self, mask, n, right=0):
        return self.rows[mask][n][right]

    def extend(self, maximum):
        self.strip.extend(maximum)
        for n in range(len(self.rows[self.core.full]), maximum + 1):
            for mask in self.core.states:
                self.rows[mask].append([mex({g for _, g in self.options(mask, n, right)})
                                       for right in range(8)])
        return self

    def options(self, mask, n, right=0, certificate=None):
        if certificate:
            start, period = certificate
            pure = lambda length, p, l, r: self.strip.lookup(length, p, l, r, start, period)
            prefix = lambda state, length, r: self.lookup(state, length, r, start, period)
            positions = representative_columns(n, start, period)
        else:
            pure, prefix = self.strip.value, self.value
            positions = range(n)
        options = []
        for v, neighbours in enumerate(self.core.adjacency):
            if not mask >> v & 1:
                continue
            rest = mask & ~self.core.closed[v]
            if not n:
                value = self.core.core_value(rest)
            elif v == 0:
                value = self.core.core_value(rest) ^ pure(n, 0, self.first, right)
            elif neighbours & 1:
                value = self.core.core_value(rest) ^ pure(n, 0, 0, right)
            else:
                value = prefix(rest, n, right)
            options.append((("core", v), value))
        for i in positions:
            playable = self.strip.word[i % self.strip.degree]
            if i == n - 1:
                playable &= ~right
            for bit in (1, 2, 4):
                if playable & bit:
                    left = (self.core.core_value(mask & ~1 if bit == self.first else mask) if i == 0
                            else prefix(mask, i, bit))
                    value = left ^ pure(n - i - 1, i + 1, bit, right)
                    options.append((("tail", bit.bit_length(), i + 1), value))
        return options

    def certificate(self, period, minimum=1):
        if period % self.strip.degree:
            raise ValueError("period must preserve the tail phase")
        rows = [self.strip.rows[n] + [g for mask in self.core.states
                                     for g in self.rows[mask][n]]
                for n in range(len(self.rows[self.core.full]))]
        return period_start(rows, period, minimum)

    def lookup(self, mask, n, right, start, period):
        index = n if n < start else start + (n - start) % period
        return self.value(mask, index, right)

    def report(self, period, start):
        if self.certificate(period, start) != start:
            raise ValueError("rooted period window does not pass")
        end = window_end(start, period)
        states = [{"mask": format(mask, "x"), "rows": self.rows[mask][:end + 1]}
                  for mask in self.core.states]
        return {"adjacency": list(self.core.adjacency), "attachment": self.first, "period": period,
                "start": start, "maximum_length": end, "states": states,
                "full_period_block": [self.value(self.core.full, n) for n in
                                      range(start, start + period)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--word", default="4646313136")
    parser.add_argument("--columns", help="solve a finite arbitrary three-row word (masks 0..7)")
    parser.add_argument("--maximum", type=int, default=80)
    parser.add_argument("--period", type=int, default=10)
    parser.add_argument("--length", type=int, default=10**12 + 1)
    parser.add_argument("--head", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--examples", action="store_true")
    args = parser.parse_args()
    if args.columns is not None:
        model = DefectiveThreeRowIntervals(map(int, args.columns))
        value = model.full_value()
        print(json.dumps({"columns": args.columns, "special_columns": len(model.special),
                          "SG": value, "interval_states": model.value.cache_info().currsize}))
        return
    if args.examples:
        report = example_report()
        raw = json.dumps(report, ensure_ascii=False, separators=(",", ":")) + "\n"
        if args.output:
            args.output.write_text(raw)
        print(json.dumps({"words": len(report["strips"]), "cores": len(report["cores"]),
                          "output": str(args.output) if args.output else None}))
        return
    if args.maximum < 1 or args.length < 0:
        parser.error("maximum must be positive and length nonnegative")
    strip = PeriodicStrip(tuple(map(int, args.word))).extend(args.maximum)
    model = strip
    core = None
    if args.head is not None:
        if args.head < 1:
            parser.error("head width must be positive")
        if strip.word[0] != 4:
            parser.error("this grid head requires row 3 in the first tail column")
        from strip3_defect_algebra import full_head
        cells, adjacency = full_head(args.head)
        core = RootedPeriodicStrip(adjacency, strip).extend(args.maximum)
        model = core
    start = model.certificate(args.period)
    if start is None:
        print(json.dumps({"status": "no finite certificate", "maximum": args.maximum,
                          "period_tested": args.period}))
        return
    if core:
        value = core.lookup(core.core.full, args.length, 0, start, args.period)
        options = core.options(core.core.full, args.length, certificate=(start, args.period))
        zero = next((move for move, g in options if g == 0), None)
        if zero:
            zero = cells[zero[1]] if zero[0] == "core" else (zero[1], args.head + zero[2])
    else:
        value = strip.lookup(args.length, 0, 0, 0, start, args.period)
        options = strip.options(args.length, certificate=(start, args.period))
        zero = next((move for move, g in options if g == 0), None)
    report = {"format": "nodekayles-clique-strip-period-v1",
              "strip": strip.report(args.period, start)}
    if core:
        report["core"] = core.report(args.period, start)
    if args.output:
        args.output.write_text(json.dumps(report, ensure_ascii=False) + "\n")
    print(json.dumps({"word": args.word, "head": args.head, "length": args.length,
                      "SG": value, "period": args.period, "start": start,
                      "representative_moves": len(options), "winning_move": zero}))


if __name__ == "__main__":
    main()
