#!/usr/bin/env python3
"""A bijective language for every reachable connected three-row component.

Column bits 1,2,4 are the top, middle, bottom rows. The index is over physical
column masks, before reflections; it is not a quotient by game equivalence.
See notes/17_connected_history_language.md for the arbitrary-width proof.
"""

import argparse
from dataclasses import dataclass
from functools import cache
from itertools import product
import json

OUTER = (0, 1, 4)
ENDS = (0, 1, 2, 4, 5)


def regular(word):
    previous = None
    for i, bit in enumerate(word):
        if bit not in OUTER:
            return False
        if bit:
            if previous is not None:
                j, side = previous
                if i - j < (2 if side == bit else 3):
                    return False
            previous = i, bit
    return True


def words(n):
    """Generate the regular language without examining rejected words."""
    def visit(prefix, side, age):
        if len(prefix) == n:
            yield prefix
            return
        for bit in OUTER:
            if not bit:
                yield from visit(prefix + (0,), side, min(4, age + 1))
            elif not side or age >= (2 if side == bit else 3):
                yield from visit(prefix + (bit,), bit, 1)
    yield from visit((), 0, 4)


def outer_mask(word):
    if not regular(word):
        raise ValueError("outer history violates the connected-language gaps")
    columns = [7] * len(word)
    for i, bit in enumerate(word):
        if bit:
            columns[i] &= ~(bit | (bit << 1) | (bit >> 1))
            if i:
                columns[i - 1] &= ~bit
            if i + 1 < len(word):
                columns[i + 1] &= ~bit
    return sum(col << (3 * i) for i, col in enumerate(columns))


def neighbourhoods(n):
    return tuple(sum(1 << (3 * d + s) for d, s in
                     ((c, r), (c, r - 1), (c, r + 1), (c - 1, r), (c + 1, r))
                     if 0 <= d < n and 0 <= s < 3)
                 for c in range(n) for r in range(3))


def components(mask, closed):
    while mask:
        reached = frontier = mask & -mask
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            new = closed[bit.bit_length() - 1] & mask & ~reached
            reached |= new
            frontier |= new
        mask ^= reached
        yield reached


def trim(mask):
    if not mask:
        return 0
    return mask >> (3 * (((mask & -mask).bit_length() - 1) // 3))


def endpoint_choices(word, right=False):
    prefix = word[::-1] if right else word
    distance = next((i for i, bit in enumerate(prefix[:3]) if bit), 3)
    return ((0,), (0, 2), (0, 1, 4, 5), ENDS)[distance]


@dataclass(frozen=True)
class Code:
    word: tuple
    left: int = 0
    right: int = 0

    def mask(self):
        if not self.word or self.left not in ENDS or self.right not in ENDS:
            raise ValueError("invalid component code")
        return outer_mask(self.word) & ~self.left & ~(self.right << (3 * (len(self.word) - 1)))

    def witness(self):
        """An explicit independent history on at most width+4 columns."""
        offset = 2 if self.left else 0
        n = len(self.word) + offset + (2 if self.right else 0)
        history = {(bit.bit_length() - 1, i + offset)
                   for i, bit in enumerate(self.word) if bit}
        for deletion, boundary, direction in (
                (self.left, offset - 1, -1),
                (self.right, offset + len(self.word), 1)):
            if deletion in (1, 4):
                row = deletion.bit_length() - 1
                history.add((row, boundary))
                history.add((2 - row, boundary + direction))
            elif deletion == 2:
                history.add((1, boundary))
            elif deletion == 5:
                history.update(((0, boundary), (2, boundary)))
        return n, tuple(sorted(history)), self.mask() << (3 * offset)


@cache
def small_catalogue(n):
    """Only widths 1/2 need a finite dictionary instead of the unique decoder."""
    if n not in (1, 2):
        raise ValueError("small catalogue is defined only for widths 1 and 2")
    closed = neighbourhoods(n)
    found = {}
    for word in words(n):
        for left, right in product(ENDS, repeat=2):
            code = Code(word, left, right)
            mask = code.mask()
            if not mask & 7 or mask.bit_length() <= 3 * (n - 1):
                continue
            if len(tuple(components(mask, closed))) == 1:
                found.setdefault(mask, code)
    return tuple(sorted(found.items()))


def decode(mask):
    """Unique normal form for width>=3; reject masks outside the family.

    This uses only the column word, not a supplied historical witness. The
    returned code reconstructs all cells and includes a legal padded witness.
    """
    if type(mask) is not int or mask <= 0 or not mask & 7:
        raise ValueError("nonempty column-trimmed mask required")
    n = (mask.bit_length() + 2) // 3
    if n < 3:
        for candidate, code in small_catalogue(n):
            if candidate == mask:
                return code
        raise ValueError("not a reachable connected component")
    columns = tuple(mask >> (3 * i) & 7 for i in range(n))
    word = [0] * n
    for i in range(1, n - 1):
        if not columns[i] & 2:
            if columns[i] not in (1, 4):
                raise ValueError("internal middle deletion must mark one outer move")
            word[i] = 5 ^ columns[i]
    for i, neighbour in ((0, 1), (n - 1, n - 2)):
        if columns[i] in (1, 4) and columns[neighbour] & 2:
            word[i] = 5 ^ columns[i]
    word = tuple(word)
    base = outer_mask(word)
    lefts, rights = endpoint_choices(word), endpoint_choices(word, True)
    left = next((d for d in lefts if (base & 7) & ~d == columns[0]), None)
    right = next((d for d in rights if (base >> (3 * (n - 1))) & ~d == columns[-1]), None)
    if left is None or right is None:
        raise ValueError("unrealizable endpoint deletion")
    code = Code(word, left, right)
    if code.mask() != mask:
        raise ValueError("columns disagree with the recovered history")
    return code


def history_count(n):
    if type(n) is not int or n < 0:
        raise ValueError("nonnegative width required")
    terms = [1, 3, 5]
    for i in range(3, n + 1):
        terms.append(sum(terms[-3:]))
    return terms[n]


def component_count(n):
    if type(n) is not int or n < 0:
        raise ValueError("nonnegative width required")
    terms = [0, 6, 24, 51, 89, 163, 301]
    for i in range(7, n + 1):
        terms.append(sum(terms[-3:]))
    return terms[n]


class LanguageIndex:
    """Perfect rank/unrank in O(n), after O(n) small automaton counts.

    Each index corresponds to one distinct mask, without a hash table or a
    stored geometry key. Widths>=3 order by outer word, left end, right end.
    """
    def __init__(self, n):
        if type(n) is not int or n < 1:
            raise ValueError("positive width required")
        self.n = n

        @cache
        def completions(pos, side, age, left_weight):
            if pos == n:
                return left_weight * (5 if not side else (1, 2, 4, 5)[age - 1])
            result = 0
            for bit, next_side, next_age, weight in self.transitions(
                    pos, side, age, left_weight):
                result += completions(pos + 1, next_side, next_age, weight)
            return result
        self.completions = completions
        self.size = (len(small_catalogue(n)) if n < 3 else completions(0, 0, 4, 0))
        assert self.size == component_count(n)

    @staticmethod
    def transitions(pos, side, age, left_weight):
        for bit in OUTER:
            if bit and side and age < (2 if bit == side else 3):
                continue
            weight = left_weight
            if pos < 3 and not weight and bit:
                weight = 1 << pos
            if pos == 2 and not weight:
                weight = 5
            yield bit, bit or side, 1 if bit else min(4, age + 1), weight

    def rank(self, mask):
        code = decode(mask)
        if len(code.word) != self.n:
            raise ValueError("mask has a different width")
        if self.n < 3:
            return next(i for i, (m, _) in enumerate(small_catalogue(self.n)) if m == mask)
        rank, side, age, weight = 0, 0, 4, 0
        for pos, chosen in enumerate(code.word):
            for bit, next_side, next_age, next_weight in self.transitions(pos, side, age, weight):
                if bit == chosen:
                    side, age, weight = next_side, next_age, next_weight
                    break
                rank += self.completions(pos + 1, next_side, next_age, next_weight)
            else:
                raise ValueError("nonregular word")
        lefts = endpoint_choices(code.word)
        rights = endpoint_choices(code.word, True)
        return rank + lefts.index(code.left) * len(rights) + rights.index(code.right)

    def unrank(self, rank):
        if type(rank) is not int or not 0 <= rank < self.size:
            raise ValueError("rank out of range")
        if self.n < 3:
            return small_catalogue(self.n)[rank][1]
        word, side, age, weight = [], 0, 4, 0
        for pos in range(self.n):
            for bit, next_side, next_age, next_weight in self.transitions(pos, side, age, weight):
                count = self.completions(pos + 1, next_side, next_age, next_weight)
                if rank < count:
                    word.append(bit)
                    side, age, weight = next_side, next_age, next_weight
                    break
                rank -= count
        word = tuple(word)
        lefts, rights = endpoint_choices(word), endpoint_choices(word, True)
        left, right = divmod(rank, len(rights))
        return Code(word, lefts[left], rights[right])


class LanguageSolver:
    """Exact move-closed SG evaluation of the entire reachable family.

    A dense table stores one byte per indexed geometry, no mask/hash key.
    This Python implementation demonstrates closure and memory accounting;
    it makes no speed claim against the specialised C solver.
    """
    def __init__(self, maximum, dense=False):
        if type(maximum) is not int or not 1 <= maximum <= 84:
            raise ValueError("supported width is 1..84 (one-byte nimbers)")
        self.maximum = maximum
        self.dense = dense
        self.closed = neighbourhoods(maximum)
        self.indices = {}
        self.tables = {}
        self.visited = 0

    def value(self, mask):
        if type(mask) is not int or mask < 0 or mask.bit_length() > 3 * self.maximum:
            raise ValueError("mask does not fit the configured board")
        result = 0
        for component in components(mask, self.closed):
            result ^= self.connected_value(trim(component))
        return result

    def connected_value(self, mask):
        n = (mask.bit_length() + 2) // 3
        if n not in self.indices:
            self.indices[n] = LanguageIndex(n)
            self.tables[n] = (bytearray([255]) * self.indices[n].size if self.dense else {})
        key = self.indices[n].rank(mask)
        table = self.tables[n]
        saved = table[key] if self.dense else table.get(key, 255)
        if saved != 255:
            return saved
        options = set()
        bits = mask
        while bits:
            bit = bits & -bits
            bits ^= bit
            options.add(self.value(mask & ~self.closed[bit.bit_length() - 1]))
        g = 0
        while g in options:
            g += 1
        assert g <= mask.bit_count() < 255
        table[key] = g
        self.visited += 1
        return g


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, default=12)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--dense", action="store_true")
    parser.add_argument("--mask", help="hexadecimal column-major connected mask")
    parser.add_argument("--rank", type=int, help="decode one perfect-index entry")
    args = parser.parse_args()
    n = args.width
    report = {"width": n, "spanning_histories": history_count(n),
              "connected_masks": component_count(n),
              "dense_sg_bytes_through_width": sum(component_count(i) for i in range(1, n + 1))}
    if args.mask is not None or args.rank is not None:
        index = LanguageIndex(n)
        code = decode(int(args.mask, 16)) if args.mask is not None else index.unrank(args.rank)
        report.update(mask=hex(code.mask()), rank=index.rank(code.mask()),
                      word="".join(str(b) for b in code.word),
                      left=code.left, right=code.right)
    if args.solve:
        solver = LanguageSolver(n, dense=args.dense)
        mask = code.mask() if args.mask is not None or args.rank is not None else (1 << (3 * n)) - 1
        report.update(grundy=solver.value(mask), evaluated_connected_masks=solver.visited)
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
