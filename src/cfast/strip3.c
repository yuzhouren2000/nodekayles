/* Exact Node Kayles on induced subgraphs of the three-row strip.
 * Column-major bits: vertex (r,c) = 3*c+r, 0 <= r < 3.
 * Up to 42 columns / 126 cells fit in unsigned __int128.
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
 * Set 0 for comparison; geometry normalisation is always enabled.
 * STRIP3_HT_PREALLOC=1 allocates the maximum table immediately instead of
 * growing from 2^18, so no doubling ever keeps two tables at once (needed
 * when the maximum already approaches the available memory).
 * STRIP3_HT_LOAD (10..95, default 70) sets the percentage of the maximum
 * capacity at which the run gives up with exit code 3.
 * The bounded table exits with code 3 when the load limit is reached.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 Mask;
enum { MAXCOLS = 42, MAXCELLS = 126 };
static Mask rowmask[3];
static Mask *keys;
static uint8_t *values;
static size_t capacity, entries, maxbits = 26, load_pct = 70;
static int reductions = 7, path_g[MAXCELLS + 1];
static uint64_t calls, moves, path_hits, cycle_hits, rotation_hits, cancelled;
static uint64_t progress_every;
static double started;

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

static uint64_t hash(Mask m) {
    uint64_t x = (uint64_t)m ^ ((uint64_t)(m >> 64) * 0x9e3779b97f4a7c15ULL);
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static void allocate(size_t bits) {
    capacity = (size_t)1 << bits;
    keys = calloc(capacity, sizeof *keys);
    values = calloc(capacity, sizeof *values);
    if (!keys || !values) {
        fprintf(stderr, "strip3: allocation failed (%zu slots)\n", capacity);
        exit(3);
    }
}

static void insert(Mask m, int g) {
    size_t slot = hash(m) & (capacity - 1);
    while (keys[slot]) slot = (slot + 1) & (capacity - 1);
    keys[slot] = m;
    values[slot] = (uint8_t)g;
}

static void store(Mask m, int g) {
    if ((entries + 1) * 100 >= capacity * load_pct) {
        if (capacity == (size_t)1 << maxbits) {
            fprintf(stderr, "strip3: memo budget reached at %zu states; incomplete, no value emitted for current board\n", entries);
            exit(3);
        }
        Mask *oldkeys = keys;
        uint8_t *oldvalues = values;
        size_t oldcap = capacity, bits = 0;
        for (size_t c = oldcap; c; c >>= 1) bits++;
        allocate(bits);
        for (size_t i = 0; i < oldcap; i++) if (oldkeys[i]) insert(oldkeys[i], oldvalues[i]);
        free(oldkeys); free(oldvalues);
    }
    insert(m, g);
    entries++;
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
            cancelled++;
        } else components[count++] = key;
    }
    int g = 0;
    for (int i = 0; i < count; i++) g ^= component_value(components[i]);
    return g;
}

static int component_value(Mask m) {
    calls++;
    if (progress_every && calls % progress_every == 0)
        fprintf(stderr, "progress calls=%llu entries=%zu %.1fs\n",
                (unsigned long long)calls, entries, seconds() - started);
    size_t slot = hash(m) & (capacity - 1);
    while (keys[slot]) {
        if (keys[slot] == m) return values[slot];
        slot = (slot + 1) & (capacity - 1);
    }
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
            if (safe) { rotation_hits++; store(m, 0); return 0; }
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
            if (endpoints || nodes == 1) { path_hits++; g = path_g[nodes]; }
            else { cycle_hits++; g = path_g[nodes - 3] == 0 ? 1 : 0; }
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
        Mask b = (Mask)1 << v;
        moves++;
        int s = value(m & ~(b | neighbours(b)));
        if (s < 128) seen[s >> 6] |= 1ULL << (s & 63);
    }
    int g = 0;
    while (g < 128 && (seen[g >> 6] >> (g & 63) & 1)) g++;
    if (g >= 128) { fprintf(stderr, "strip3: impossible mex overflow\n"); exit(1); }
    store(m, g);
    return g;
}

static Mask full_board(int n) { return ((Mask)1 << (3 * n)) - 1; }

static void stats(int n, int g, double t) {
    fprintf(stderr, "[3x%d] G=%d time=%.3f memo=%zu calls=%llu moves=%llu path=%llu cycle=%llu rotation=%llu cancelled=%llu\n",
        n, g, seconds() - t, entries, (unsigned long long)calls, (unsigned long long)moves,
        (unsigned long long)path_hits, (unsigned long long)cycle_hits,
        (unsigned long long)rotation_hits, (unsigned long long)cancelled);
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
        if (!*env || *end || flags < 0 || flags > 7) { fprintf(stderr, "invalid STRIP3_REDUCTIONS (0..7)\n"); exit(2); }
        reductions = (int)flags;
    }
    env = getenv("STRIP3_HT_PREALLOC");
    int prealloc = env && strcmp(env, "0") != 0;
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
            printf("%d\n", value(m));
            fflush(stdout);   /* keep completed lines on timeout kills */
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
    fprintf(stderr, "done %.3fs memo=%zu capacity=%zu bytes_per_slot=17 reductions=%d\n", seconds() - started, entries, capacity, reductions);
    free(keys); free(values);
    return 0;
usage:
    fprintf(stderr, "usage: strip3 strip N | survey N | board N | masks (N=1..42)\n");
    return 2;
}
