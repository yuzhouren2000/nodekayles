/* Exact Node Kayles on induced subgraphs of the three-row strip.
 * Column-major bits: vertex (r,c) = 3*c+r, 0 <= r < 3.
 * Up to 40 columns / 120 cells fit in the low 120 bits of unsigned __int128;
 * the top 8 bits of each hash-table word hold the nimber (0..255), so one
 * 16-byte word is both the key and the value: a memo hit touches a single
 * cache line and the table needs no separate value array.
 *
 * Structural reductions (each preserves the exact nimber):
 * - any connected path, including bent embeddings, uses the path recurrence;
 * - connected degree-2 cycles have one option type, a path of order k-3;
 * - equal canonical components cancel before recursive evaluation;
 * - a safe 180-degree involution gives nimber zero.
 *
 * Usage: strip3 strip N | survey N | masks | board N
 * survey emits JSON lines with G(3xn), A(n) (one end-middle deleted), D(n)
 * (both end-middles deleted), and every first-move option. masks reads decimal
 * N and hexadecimal column-major mask, one position per line.
 * Env: STRIP3_HT_BITS (max, default 26), STRIP3_REDUCTIONS is a bit set:
 * 1=paths/cycles, 2=safe rotation, 4=equal component cancellation (default 7).
 * 8=certified full single-thick-block combs (widths 2..5),
 * 16=certified full double-thick-block combs (widths 2,3, odd bridge).
 * Set 0 for comparison; geometry normalisation is always enabled.
 * STRIP3_HT_PREALLOC=1 allocates the maximum table immediately instead of
 * growing from 2^18, so no doubling ever keeps two tables at once (needed
 * when the maximum already approaches the available memory).
 * STRIP3_HT_LOAD (10..95, default 70) sets the percentage of the maximum
 * capacity at which the run gives up with exit code 3.
 * The bounded table exits with code 3 when the load limit is reached.
 */

#include <pthread.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "strip3_comb_tables.h"

typedef unsigned __int128 Mask;
enum { MAXCOLS = 42, MAXCELLS = 126 };
_Static_assert(MAXCOLS == STRIP3_COMB_MAXCOLS, "regenerate bounded comb tables");
static const Mask KEYMASK = (((Mask)1 << 120) - 1);
static Mask rowmask[3];
static Mask *slots;                 /* packed key (low 120 bits) + value (top 8) */
static size_t capacity, maxbits = 26, load_pct = 70;
static _Atomic size_t entries;      /* global count, atomic: many threads insert */
static int reductions = 7, path_g[MAXCELLS + 1];
static uint64_t progress_every;
static double started;

/* Sharded memo: the table is split into lane_count lanes of lane_size slots.
 * A key stays inside its lane (lane = low hash bits, local probe from the
 * remaining bits); each lane has a spinlock, so lookups and inserts by
 * different threads only serialise when they touch the same lane. */
static size_t lane_count = 1, lane_size = 1, lane_shift = 0;
static unsigned char *lane_lock;
static int nthreads = 1;

/* Keys of 121..126 bits (n = 41, 42; exercised by the sparse regression
 * tests, never by production runs up to n = 36) do not fit in the packed
 * 16-byte word, so they live in a separate small table with a value array. */
static Mask *wide_keys;
static uint8_t *wide_values;
static size_t wide_cap, wide_entries;

/* per-thread statistics, summed at exit */
struct tstat {
    uint64_t calls, moves, path_hits, cycle_hits, rotation_hits, cancelled;
    uint64_t single_comb, chain_comb, chain_formula;
};
static __thread struct tstat tst;
enum { TMAX = 64 };
static struct tstat stats_out[TMAX];   /* copied by each thread before it exits */
static int nregistered;

static double seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static int lowbit(Mask m) {
    uint64_t low = (uint64_t)m;
    return low ? __builtin_ctzll(low) : 64 + __builtin_ctzll((uint64_t)(m >> 64));
}

static int order(Mask m) {
    return __builtin_popcountll((uint64_t)m) + __builtin_popcountll((uint64_t)(m >> 64));
}

static Mask neighbours(Mask m) {
    return (m << 3) | (m >> 3)
        | ((m & (rowmask[0] | rowmask[1])) << 1)
        | ((m & (rowmask[1] | rowmask[2])) >> 1);
}

static uint64_t reverse64(uint64_t x) {
    x = ((x >> 1) & 0x5555555555555555ULL) | ((x & 0x5555555555555555ULL) << 1);
    x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
    x = ((x >> 4) & 0x0f0f0f0f0f0f0f0fULL) | ((x & 0x0f0f0f0f0f0f0f0fULL) << 4);
    return __builtin_bswap64(x);
}

static Mask reverse_bits(Mask m, int length) {
    Mask rev = ((Mask)reverse64((uint64_t)m) << 64) | reverse64((uint64_t)(m >> 64));
    return rev >> (128 - length);
}

static Mask flip_rows(Mask m) {
    return ((m & rowmask[0]) << 2) | (m & rowmask[1]) | ((m & rowmask[2]) >> 2);
}

static int width(Mask m) {
    uint64_t high = (uint64_t)(m >> 64);
    int last = high ? 127 - __builtin_clzll(high) : 63 - __builtin_clzll((uint64_t)m);
    return last / 3 + 1;
}

/* Trim empty columns and rows, then minimise the rectangle symmetries.
 * In height 2, reversing all 3*n bits reflects into rows 1,2; shift down.
 * A 3x3 box also permits transposition, so include all eight transforms.
 */
static Mask canonical(Mask m) {
    m >>= (lowbit(m) / 3) * 3;
    while (!(m & rowmask[0])) m >>= 1;
    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask best = m;
    for (int trans = 0; trans < ((h == 3 && n == 3) ? 2 : 1); trans++) {
        Mask base = m;
        if (trans) {
            base = 0;
            for (Mask x = m; x; x &= x - 1) {
                int v = lowbit(x);
                base |= (Mask)1 << (3 * (v % 3) + v / 3);
            }
        }
        Mask vertical = flip_rows(base) >> (3 - h);
        Mask both = reverse_bits(base, 3 * n) >> (3 - h);
        Mask horizontal = flip_rows(both) >> (3 - h);
        if (base < best) best = base;
        if (vertical < best) best = vertical;
        if (both < best) best = both;
        if (horizontal < best) best = horizontal;
    }
    return best;
}

static int column(Mask m, int c) { return (int)(m >> (3 * c)) & 7; }

static size_t comb_pair_index(int budget, int left, int right) {
    return (size_t)(left * (budget + 1) - left * (left - 1) / 2 + right);
}

/* Match the ENTIRE component, never a subword or a damaged thick block.
 * Full columns locate the blocks; caps must be 6 and every remaining column
 * must be the exact phase-zero 4,6 comb seen from its cap. Vertical reflection
 * is tried by the caller. A reversed 3,2 chain swaps the outer lengths.
 */
static int oriented_comb(Mask m, int n, Mask full, int total, int *g, int *formula) {
    int first = lowbit(full) / 3, after = first;
    while (after < n && column(m, after) == 7) after++;
    int w = after - first, left = first - 1;
    if (left < 0 || after >= n || column(m, first - 1) != 6 || column(m, after) != 6)
        return 0;
    for (int i = 0; i < left; i++)
        if (column(m, left - 1 - i) != (i % 2 ? 6 : 4)) return 0;
    if ((reductions & 8) && total == w && w >= 2 && w <= 5) {
        int right = n - after - 1;
        for (int i = 0; i < right; i++)
            if (column(m, after + 1 + i) != (i % 2 ? 6 : 4)) return 0;
        int budget = MAXCOLS - w - 2;
        *g = strip3_single_values[strip3_single_offsets[w - 2] + comb_pair_index(budget, left, right)];
        return 8;
    }
    if (!(reductions & 16) || total != 5 || (w != 2 && w != 3)) return 0;
    int second = after + 1;
    while (second < n && column(m, second) != 7) second++;
    int bridge = second - after - 2, end = second;
    if (second >= n || bridge < 1 || !(bridge % 2) || column(m, second - 1) != 6) return 0;
    for (int i = 0; i < bridge; i++)
        if (column(m, after + 1 + i) != (i % 2 ? 6 : 4)) return 0;
    while (end < n && column(m, end) == 7) end++;
    if (end - second != 5 - w || end >= n || column(m, end) != 6) return 0;
    int right = n - end - 1;
    for (int i = 0; i < right; i++)
        if (column(m, end + 1 + i) != (i % 2 ? 6 : 4)) return 0;
    if (w == 3) { int temp = left; left = right; right = temp; }
    int budget = MAXCOLS - 9 - bridge;
    *g = strip3_chain_values[strip3_chain_offsets[(bridge - 1) / 2] + comb_pair_index(budget, left, right)];
    *formula = left >= 5 && right >= 5 && bridge >= 11;
    return 16;
}

static int certified_comb(Mask m, int n, int *g, int *formula) {
    Mask full = m & (m >> 1) & (m >> 2) & rowmask[0];
    if (!full) return 0;
    int total = order(full);
    if (!(((reductions & 8) && total >= 2 && total <= 5) || ((reductions & 16) && total == 5)))
        return 0;
    int kind = oriented_comb(m, n, full, total, g, formula);
    return kind ? kind : oriented_comb(flip_rows(m), n, full, total, g, formula);
}

static uint64_t hash(Mask m) {
    uint64_t x = (uint64_t)m ^ ((uint64_t)(m >> 64) * 0x9e3779b97f4a7c15ULL);
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static inline void lane_acquire(size_t lane) {
    if (lane_count == 1) return;
    unsigned char *l = &lane_lock[lane];
    while (__atomic_test_and_set(l, __ATOMIC_ACQUIRE)) {
        while (__atomic_load_n(l, __ATOMIC_RELAXED)) __asm__ volatile("pause");
    }
}
static inline void lane_release(size_t lane) {
    if (lane_count == 1) return;
    __atomic_clear(&lane_lock[lane], __ATOMIC_RELEASE);
}

static void allocate(size_t bits) {
    capacity = (size_t)1 << bits;
    while (lane_count > 1 && capacity / lane_count < 1024) lane_count >>= 1;
    lane_size = capacity / lane_count;
    lane_shift = 0;
    for (size_t c = lane_count; c > 1; c >>= 1) lane_shift++;
    slots = calloc(capacity, sizeof *slots);
    if (!slots) {
        fprintf(stderr, "strip3: allocation failed (%zu slots)\n", capacity);
        exit(3);
    }
    if (lane_count > 1) {
        lane_lock = calloc(lane_count, 1);
        if (!lane_lock) { fprintf(stderr, "strip3: lock allocation failed\n"); exit(3); }
    }
    if (!wide_keys) {
        wide_cap = (size_t)1 << 22;
        wide_keys = calloc(wide_cap, sizeof *wide_keys);
        wide_values = calloc(wide_cap, sizeof *wide_values);
        if (!wide_keys || !wide_values) { fprintf(stderr, "strip3: wide table allocation failed\n"); exit(3); }
        wide_entries = 0;
    }
}

/* wide-key memo: same lane locks, separate table */
static int wide_get(Mask m) {
    size_t hh = hash(m);
    size_t lane = lane_count == 1 ? 0 : hh & (lane_count - 1);
    lane_acquire(lane);
    size_t slot = hh & (wide_cap - 1);
    for (;;) {
        if (!wide_keys[slot]) { lane_release(lane); return -1; }
        if (wide_keys[slot] == m) { int g = wide_values[slot]; lane_release(lane); return g; }
        slot = (slot + 1) & (wide_cap - 1);
    }
}
static void wide_put(Mask m, int g) {
    size_t hh = hash(m);
    size_t lane = lane_count == 1 ? 0 : hh & (lane_count - 1);
    lane_acquire(lane);
    size_t slot = hh & (wide_cap - 1);
    while (wide_keys[slot] && wide_keys[slot] != m) slot = (slot + 1) & (wide_cap - 1);
    if (!wide_keys[slot]) {
        if ((wide_entries + 1) * 10 >= wide_cap * 7) { lane_release(lane); fprintf(stderr, "strip3: wide memo full\n"); exit(3); }
        wide_keys[slot] = m; wide_values[slot] = (uint8_t)g;
        __atomic_fetch_add(&wide_entries, 1, __ATOMIC_RELAXED);
    }
    lane_release(lane);
}

static void insert(Mask m, int g) {
    size_t h = hash(m);
    size_t lane = lane_count == 1 ? 0 : h & (lane_count - 1);
    size_t base = lane * lane_size;
    size_t local = lane_count == 1 ? (h & (lane_size - 1))
                                   : ((h >> lane_shift) & (lane_size - 1));
    Mask entry = m | ((Mask)(uint8_t)g << 120);
    lane_acquire(lane);
    size_t step = 0;
    for (;;) {
        Mask cur = slots[base + local];
        if (!(cur & KEYMASK)) { slots[base + local] = entry; __atomic_fetch_add(&entries, 1, __ATOMIC_RELAXED); break; }
        if ((cur & KEYMASK) == m) break;          /* raced insert: value is identical */
        local = (local + 1) & (lane_size - 1);
        if (++step > lane_size) { lane_release(lane); fprintf(stderr, "strip3: lane overflow\n"); exit(3); }
    }
    lane_release(lane);
}

static void store(Mask m, int g) {
    size_t count = __atomic_load_n(&entries, __ATOMIC_RELAXED);
    if ((count + 1) * 100 >= capacity * load_pct) {
        if (capacity == (size_t)1 << maxbits) {
            fprintf(stderr, "strip3: memo budget reached at %zu states; incomplete, no value emitted for current board\n", count);
            exit(3);
        }
        if (nthreads > 1) { fprintf(stderr, "strip3: growth forbidden with threads>1 (use STRIP3_HT_PREALLOC=1)\n"); exit(3); }
        Mask *oldslots = slots;
        size_t oldcap = capacity, bits = 0;
        for (size_t c = oldcap; c; c >>= 1) bits++;
        allocate(bits);
        __atomic_store_n(&entries, (size_t)0, __ATOMIC_RELAXED);
        for (size_t i = 0; i < oldcap; i++) if (oldslots[i] & KEYMASK) {
            Mask m2 = oldslots[i] & KEYMASK;
            int g2 = (int)(oldslots[i] >> 120);
            insert(m2, g2);
        }
        free(oldslots);
    }
    if (m >> 120) { wide_put(m, g); return; }
    insert(m, g);
}

static int value(Mask m);
static int component_value(Mask m);

static int value(Mask m) {
    Mask components[MAXCELLS];
    int count = 0;
    while (m) {
        Mask comp = (Mask)1 << lowbit(m), frontier = comp;
        while (frontier) {
            frontier = neighbours(frontier) & m & ~comp;
            comp |= frontier;
        }
        m &= ~comp;
        Mask key = canonical(comp);
        int duplicate = -1;
        if (reductions & 4)
            for (int i = 0; i < count; i++) if (components[i] == key) { duplicate = i; break; }
        if (duplicate >= 0) {
            components[duplicate] = components[--count];
            tst.cancelled++;
        } else components[count++] = key;
    }
    int g = 0;
    for (int i = 0; i < count; i++) g ^= component_value(components[i]);
    return g;
}

static int component_value(Mask m) {
    tst.calls++;
    if (progress_every && tst.calls % progress_every == 0)
        fprintf(stderr, "progress calls=%llu entries=%zu %.1fs\n",
                (unsigned long long)tst.calls, __atomic_load_n(&entries, __ATOMIC_RELAXED),
                seconds() - started);
    if (m >> 120) {
        int g = wide_get(m);
        if (g >= 0) return g;
    }
    size_t hh = hash(m);
    size_t lane = lane_count == 1 ? 0 : hh & (lane_count - 1);
    size_t base = lane * lane_size;
    size_t local = lane_count == 1 ? (hh & (lane_size - 1))
                                   : ((hh >> lane_shift) & (lane_size - 1));
    lane_acquire(lane);
    for (;;) {
        Mask entry = slots[base + local];
        if (!(entry & KEYMASK)) break;
        if ((entry & KEYMASK) == m) { int g = (int)(entry >> 120); lane_release(lane); return g; }
        local = (local + 1) & (lane_size - 1);
    }
    lane_release(lane);
    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask vertical = flip_rows(m) >> (3 - h);
    Mask both = reverse_bits(m, 3 * n) >> (3 - h);
    Mask horizontal = flip_rows(both) >> (3 - h);
    if (reductions) {
        if ((reductions & 2) && m == both) {
            int safe = 1;
            for (Mask x = m; x && safe; x &= x - 1) {
                int v = lowbit(x), partner = 3 * (n - 1 - v / 3) + h - 1 - v % 3;
                Mask b = (Mask)1 << v;
                if (partner == v || (neighbours(b) & ((Mask)1 << partner))) safe = 0;
            }
            if (safe) { tst.rotation_hits++; store(m, 0); return 0; }
        }
        int maxdegree = 0, endpoints = 0, nodes = order(m);
        for (Mask x = (reductions & 1) ? m : 0; x; x &= x - 1) {
            int degree = order(neighbours((Mask)1 << lowbit(x)) & m);
            if (degree > maxdegree) maxdegree = degree;
            if (degree < 2) endpoints++;
            if (maxdegree > 2) break;
        }
        if ((reductions & 1) && maxdegree <= 2) {
            int g;
            if (endpoints || nodes == 1) { tst.path_hits++; g = path_g[nodes]; }
            else { tst.cycle_hits++; g = path_g[nodes - 3] == 0 ? 1 : 0; }
            store(m, g);
            return g;
        }

    }
    if (reductions & 24) {
        int g, formula = 0, kind = certified_comb(m, n, &g, &formula);
        if (kind) {
            if (kind == 8) tst.single_comb++;
            else { tst.chain_comb++; tst.chain_formula += formula; }
            store(m, g);
            return g;
        }
    }
    uint64_t seen[2] = {0, 0};
    for (Mask x = m; x; x &= x - 1) {
        int v = lowbit(x), col = v / 3, row = v % 3;
        /* Move orbits under preserving rectangle symmetries. */
        if (m == vertical && 3 * col + h - 1 - row < v) continue;
        if (m == horizontal && 3 * (n - 1 - col) + row < v) continue;
        if (m == both && 3 * (n - 1 - col) + h - 1 - row < v) continue;
        tst.moves++;
        Mask kill = (Mask)1 << v;
        if (v >= 3) kill |= (Mask)1 << (v - 3);
        if (v + 3 < 3 * n) kill |= (Mask)1 << (v + 3);
        if (row > 0) kill |= (Mask)1 << (v - 1);
        if (row < 2) kill |= (Mask)1 << (v + 1);
        int s = value(m & ~kill);
        if (s < 128) seen[s >> 6] |= 1ULL << (s & 63);
    }
    int g = 0;
    while (g < 128 && (seen[g >> 6] >> (g & 63) & 1)) g++;
    if (g >= 128) { fprintf(stderr, "strip3: impossible mex overflow\n"); exit(1); }
    store(m, g);
    return g;
}

static Mask *job_mask;
static int *job_result;
static unsigned char *job_done;
static int job_count, job_cap;
static _Atomic int job_next;
static _Atomic int job_printed;
static unsigned char print_lock;

static void *worker_main(void *arg) {
    int worker = (int)(intptr_t)arg;
    for (;;) {
        int i = __atomic_fetch_add(&job_next, 1, __ATOMIC_RELAXED);
        if (i >= job_count) break;
        job_result[i] = value(job_mask[i]);
        __atomic_store_n(&job_done[i], 1, __ATOMIC_RELEASE);
        /* flush the longest completed prefix, in order */
        while (__atomic_test_and_set(&print_lock, __ATOMIC_ACQUIRE))
            __asm__ volatile("pause");
        for (;;) {
            int k = __atomic_load_n(&job_printed, __ATOMIC_RELAXED);
            if (k >= job_count || !__atomic_load_n(&job_done[k], __ATOMIC_ACQUIRE)) break;
            printf("%d\n", job_result[k]);
            fflush(stdout);
            __atomic_store_n(&job_printed, k + 1, __ATOMIC_RELAXED);
        }
        __atomic_clear(&print_lock, __ATOMIC_RELEASE);
    }
    stats_out[worker] = tst;   /* copy before the thread's TLS disappears */
    return NULL;
}

static Mask full_board(int n) { return ((Mask)1 << (3 * n)) - 1; }

static void stats(int n, int g, double t) {
    fprintf(stderr, "[3x%d] G=%d time=%.3f memo=%zu calls=%llu moves=%llu path=%llu cycle=%llu rotation=%llu cancelled=%llu single_comb=%llu chain_comb=%llu chain_formula=%llu\n",
        n, g, seconds() - t, __atomic_load_n(&entries, __ATOMIC_RELAXED),
        (unsigned long long)tst.calls, (unsigned long long)tst.moves,
        (unsigned long long)tst.path_hits, (unsigned long long)tst.cycle_hits,
        (unsigned long long)tst.rotation_hits, (unsigned long long)tst.cancelled,
        (unsigned long long)tst.single_comb, (unsigned long long)tst.chain_comb, (unsigned long long)tst.chain_formula);
}

static void initialise(void) {
    for (int c = 0; c < MAXCOLS; c++)
        for (int r = 0; r < 3; r++) rowmask[r] |= (Mask)1 << (3 * c + r);
    for (int size = 1; size <= MAXCELLS; size++) {
        uint64_t seen[2] = {0, 0};
        for (int i = 0; i < size; i++) {
            int left = i > 0 ? i - 1 : 0, right = size - i > 1 ? size - i - 2 : 0;
            int s = path_g[left] ^ path_g[right];
            if (s < 128) seen[s >> 6] |= 1ULL << (s & 63);
        }
        while (seen[path_g[size] >> 6] >> (path_g[size] & 63) & 1) path_g[size]++;
    }
    const char *env = getenv("STRIP3_HT_BITS");
    if (env) {
        char *end;
        long bits = strtol(env, &end, 10);
        if (!*env || *end || bits < 8 || bits > 30) { fprintf(stderr, "invalid STRIP3_HT_BITS (8..30)\n"); exit(2); }
        maxbits = (size_t)bits;
    }
    env = getenv("STRIP3_REDUCTIONS");
    if (env) {
        char *end;
        long flags = strtol(env, &end, 10);
        if (!*env || *end || flags < 0 || flags > 31) { fprintf(stderr, "invalid STRIP3_REDUCTIONS (0..31)\n"); exit(2); }
        reductions = (int)flags;
    }
    env = getenv("STRIP3_HT_PREALLOC");
    int prealloc = env && strcmp(env, "0") != 0;
    env = getenv("STRIP3_THREADS");
    if (env) {
        char *end;
        long n = strtol(env, &end, 10);
        if (!*env || *end || n < 1 || n > TMAX) { fprintf(stderr, "invalid STRIP3_THREADS (1..%d)\n", TMAX); exit(2); }
        nthreads = (int)n;
    }
    if (nthreads > 1 && !prealloc) {
        fprintf(stderr, "strip3: STRIP3_THREADS>1 requires STRIP3_HT_PREALLOC=1\n");
        exit(2);
    }
    if (nthreads > 1) lane_count = 1 << 16;
    env = getenv("STRIP3_HT_LOAD");
    if (env) {
        char *end;
        long pct = strtol(env, &end, 10);
        if (!*env || *end || pct < 10 || pct > 95) { fprintf(stderr, "invalid STRIP3_HT_LOAD (10..95)\n"); exit(2); }
        load_pct = (size_t)pct;
    }
    env = getenv("STRIP3_PROGRESS");
    if (env) {
        char *end;
        long every = strtol(env, &end, 10);
        if (!*env || *end || every < 0) { fprintf(stderr, "invalid STRIP3_PROGRESS (>=0)\n"); exit(2); }
        progress_every = (uint64_t)every;
    }
    allocate(prealloc || maxbits < 18 ? maxbits : 18);
    started = seconds();
}

int main(int argc, char **argv) {
    if (argc < 2) goto usage;
    if (!strcmp(argv[1], "masks") && argc == 2) {
        initialise();
        int n;
        char hex[40], line[128], extra;
        int njobs = 0;
        while (fgets(line, sizeof line, stdin)) {
            if (sscanf(line, "%d %39s %c", &n, hex, &extra) != 2) goto usage;
            if (n < 1 || n > MAXCOLS || strlen(hex) > 32) goto usage;
            Mask m = 0;
            for (const char *p = hex; *p; p++) {
                int digit = *p >= '0' && *p <= '9' ? *p - '0' :
                    *p >= 'a' && *p <= 'f' ? *p - 'a' + 10 :
                    *p >= 'A' && *p <= 'F' ? *p - 'A' + 10 : -1;
                if (digit < 0) goto usage;
                m = (m << 4) | digit;
            }
            if (m & ~full_board(n)) goto usage;
            if (njobs == job_cap) {
                int cap = job_cap ? job_cap * 2 : 1024;
                Mask *nm = realloc(job_mask, (size_t)cap * sizeof *nm);
                int *nr = realloc(job_result, (size_t)cap * sizeof *nr);
                unsigned char *nd = realloc(job_done, (size_t)cap * sizeof *nd);
                if (!nm || !nr || !nd) { fprintf(stderr, "strip3: input too large\n"); return 2; }
                job_mask = nm; job_result = nr; job_done = nd;
                job_cap = cap;
                memset(job_done + njobs, 0, (size_t)(cap - njobs));
            }
            job_mask[njobs++] = m;
        }
        if (nthreads <= 1 || njobs <= 1) {
            nregistered = 1;
            for (int i = 0; i < njobs; i++) {
                printf("%d\n", value(job_mask[i]));
                fflush(stdout);
            }
            stats_out[0] = tst;
        } else {
            int nw = njobs < nthreads ? njobs : nthreads;
            nregistered = nw;
            job_count = njobs;
            __atomic_store_n(&job_next, 0, __ATOMIC_RELAXED);
            __atomic_store_n(&job_printed, 0, __ATOMIC_RELAXED);
            pthread_t th[TMAX];
            for (int w = 1; w < nw; w++)
                pthread_create(&th[w], NULL, worker_main, (void *)(intptr_t)w);
            worker_main((void *)(intptr_t)0);
            for (int w = 1; w < nw; w++) pthread_join(th[w], NULL);
            /* anything the print gate has not flushed (cannot happen, but be safe) */
            for (int i = 0; i < njobs; i++) if (!job_done[i]) printf("%d\n", job_result[i]);
        }
    } else if (argc == 3 && (!strcmp(argv[1], "strip") || !strcmp(argv[1], "survey") || !strcmp(argv[1], "board"))) {
        char *end;
        long limit = strtol(argv[2], &end, 10);
        if (!*argv[2] || *end || limit < 1 || limit > MAXCOLS) goto usage;
        initialise();
        int first = !strcmp(argv[1], "board") ? (int)limit : 1;
        for (int n = first; n <= limit; n++) {
            double t = seconds();
            Mask full = full_board(n);
            int g = value(full);
            if (!strcmp(argv[1], "survey")) {
                Mask end_gap = full & ~((Mask)1 << (3 * (n - 1) + 1));
                int a = value(end_gap);
                int d = value(end_gap & ~((Mask)1 << 1));
                int options[MAXCELLS];
                for (int v = 0; v < 3 * n; v++) {
                    Mask b = (Mask)1 << v;
                    options[v] = value(full & ~(b | neighbours(b)));
                }
                printf("{\"n\":%d,\"grundy\":%d,\"end_gap\":%d,\"both_gaps\":%d,\"options\":[", n, g, a, d);
                for (int v = 0; v < 3 * n; v++) printf("%s%d", v ? "," : "", options[v]);
                printf("]}\n");
            } else printf("3\t%d\t%d\n", n, g);
            fflush(stdout);
            stats(n, g, t);
        }
    } else goto usage;
    {
        uint64_t tc = 0, tm = 0, tp = 0, tcy = 0, tr = 0, tcn = 0, tsc = 0, tcc = 0, tcf = 0;
        int reg = nregistered ? nregistered : 0;
        if (!reg) { stats_out[0] = tst; reg = 1; }   /* survey/strip/board modes */
        for (int i = 0; i < reg && i < TMAX; i++) {
            tc += stats_out[i].calls; tm += stats_out[i].moves;
            tp += stats_out[i].path_hits; tcy += stats_out[i].cycle_hits;
            tr += stats_out[i].rotation_hits; tcn += stats_out[i].cancelled;
            tsc += stats_out[i].single_comb; tcc += stats_out[i].chain_comb; tcf += stats_out[i].chain_formula;
        }
        fprintf(stderr, "done %.3fs memo=%zu capacity=%zu bytes_per_slot=16 reductions=%d threads=%d lanes=%zu\n",
                seconds() - started, __atomic_load_n(&entries, __ATOMIC_RELAXED), capacity, reductions, nthreads, lane_count);
        fprintf(stderr, "stats calls=%llu moves=%llu path=%llu cycle=%llu rotation=%llu cancelled=%llu single_comb=%llu chain_comb=%llu chain_formula=%llu\n",
                (unsigned long long)tc, (unsigned long long)tm, (unsigned long long)tp,
                (unsigned long long)tcy, (unsigned long long)tr, (unsigned long long)tcn,
                (unsigned long long)tsc, (unsigned long long)tcc, (unsigned long long)tcf);
    }
    free(slots);
    return 0;
usage:
    fprintf(stderr, "usage: strip3 strip N | survey N | board N | masks (N=1..42)\n");
    return 2;
}
