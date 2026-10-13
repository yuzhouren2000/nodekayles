/* Exact Node Kayles restricted to the entire reachable-component family of
 * the full three-row board. Perfect indexing: one byte per geometry, no keys.
 * See notes/17_connected_history_language.md. This is single-threaded.
 *
 * Usage: strip3_language board N | strip N | masks | ranks
 * STRIP3_LANGUAGE_BYTES caps allocated SG bytes (default 1 GiB).
 * STRIP3_LANGUAGE_REDUCTIONS uses the original 1/2/4 reduction flags.
 * masks rejects connected geometries outside the proved reachable family.
 *
 * Reuse only the original geometry primitives and path table. Its hash-table
 * initialiser/CLI are never called by this entry point.
 */
#define main original_strip3_main
#include "strip3.c"
#undef main
#include <errno.h>

static uint8_t *dense[MAXCOLS + 1];
static uint64_t ways[MAXCOLS + 1][MAXCOLS + 1][3][5][6];
static uint64_t state_counts[MAXCOLS + 1];
static uint64_t allocated_bytes, byte_budget = (uint64_t)1 << 30;
static uint64_t language_entries;
static int language_flags = 7;
static const int end_choices[4][5] = {{0}, {0,2}, {0,1,4,5}, {0,1,2,4,5}};
static const int end_sizes[4] = {1,2,4,5};
static const unsigned small1[6] = {1,2,3,4,6,7};
static const unsigned small2[24] = {
    011,013,022,023,026,027,031,032,033,036,037,044,
    046,057,062,063,064,066,067,072,073,075,076,077
};

static void invalid_mask(Mask m) {
    fprintf(stderr, "strip3_language: unreachable connected mask %llx%016llx\n",
            (unsigned long long)(m >> 64), (unsigned long long)m);
    exit(2);
}

/* Prefix automaton: side 0/1/2; age since previous outer move, capped at 4.
 * At length 3 the left boundary weight is frozen; at the end multiply by
 * the right boundary weight. Counts are <=2^64 for the supported widths.
 */
static uint64_t completion(int n, int pos, int side, int age, int weight) {
    if (pos == n) return (uint64_t)weight * (side ? end_sizes[age - 1] : 5);
    uint64_t *saved = &ways[n][pos][side][age][weight];
    if (*saved) return *saved;
    uint64_t total = 0;
    for (int bit = 0; bit < 3; bit++) {
        if (bit && side && age < (side == bit ? 2 : 3)) continue;
        int w = weight;
        if (pos < 3 && !w && bit) w = 1 << pos;
        if (pos == 2 && !w) w = 5;
        total += completion(n, pos + 1, bit ? bit : side,
                            bit ? 1 : (age < 4 ? age + 1 : 4), w);
    }
    *saved = total;
    return total;
}

static uint64_t language_rank(Mask m) {
    int n = width(m);
    if (n < 3) {
        const unsigned *catalogue = n == 1 ? small1 : small2;
        int count = n == 1 ? 6 : 24;
        for (int i = 0; i < count; i++) if (m == catalogue[i]) return (uint64_t)i;
        invalid_mask(m);
    }
    int col[MAXCOLS], word[MAXCOLS] = {0};
    for (int i = 0; i < n; i++) col[i] = (int)((m >> (3*i)) & 7);
    for (int i = 1; i < n - 1; i++) if (!(col[i] & 2)) {
        if (col[i] != 1 && col[i] != 4) invalid_mask(m);
        word[i] = col[i] == 1 ? 2 : 1;
    }
    if ((col[0] == 1 || col[0] == 4) && (col[1] & 2)) word[0] = col[0] == 1 ? 2 : 1;
    if ((col[n-1] == 1 || col[n-1] == 4) && (col[n-2] & 2)) word[n-1] = col[n-1] == 1 ? 2 : 1;
    Mask base = full_board(n);
    for (int i = 0; i < n; i++) if (word[i]) {
        Mask bit = (Mask)1 << (3*i + (word[i] == 1 ? 0 : 2));
        base &= ~(bit | neighbours(bit));
    }
    if (((base ^ m) & ~((Mask)7 | ((Mask)7 << (3*(n-1))))) || (m & ~base)) invalid_mask(m);
    int dl = 3, dr = 3;
    for (int i = 0; i < 3; i++) if (word[i]) { dl = i; break; }
    for (int i = 0; i < 3; i++) if (word[n-1-i]) { dr = i; break; }
    int il = -1, ir = -1;
    for (int i = 0; i < end_sizes[dl]; i++)
        if (((int)(base & 7) & ~end_choices[dl][i]) == col[0]) il = i;
    for (int i = 0; i < end_sizes[dr]; i++)
        if (((int)((base >> (3*(n-1))) & 7) & ~end_choices[dr][i]) == col[n-1]) ir = i;
    if (il < 0 || ir < 0) invalid_mask(m);
    uint64_t rank = 0;
    int side = 0, age = 4, weight = 0;
    for (int pos = 0; pos < n; pos++) {
        int found = 0;
        for (int bit = 0; bit < 3; bit++) {
            if (bit && side && age < (side == bit ? 2 : 3)) continue;
            int ns = bit ? bit : side, na = bit ? 1 : (age < 4 ? age + 1 : 4), nw = weight;
            if (pos < 3 && !nw && bit) nw = 1 << pos;
            if (pos == 2 && !nw) nw = 5;
            if (bit == word[pos]) { side = ns; age = na; weight = nw; found = 1; break; }
            rank += completion(n, pos + 1, ns, na, nw);
        }
        if (!found) invalid_mask(m);
    }
    rank += (uint64_t)il * end_sizes[dr] + ir;
    if (rank >= state_counts[n]) invalid_mask(m);
    return rank;
}

/* Preserve the physical three rows. Trimming rows or transposing arbitrary
 * graph embeddings can leave the reachable language and is unnecessary.
 */
static Mask language_canonical(Mask m) {
    m >>= 3 * (lowbit(m) / 3);
    int n = width(m);
    Mask best = m, vertical = flip_rows(m), both = reverse_bits(m, 3*n);
    Mask horizontal = flip_rows(both);
    if (vertical < best) best = vertical;
    if (both < best) best = both;
    if (horizontal < best) best = horizontal;
    return best;
}

static int language_value(Mask m);

static void save_language(int n, uint64_t rank, int g) {
    dense[n][rank] = (uint8_t)(g + 1);
    language_entries++;
}

static int language_component(Mask m) {
    tst.calls++;
    int n = width(m);
    uint64_t rank = language_rank(m);
    if (!dense[n]) {
        if (state_counts[n] > byte_budget - allocated_bytes) {
            fprintf(stderr, "strip3_language: byte budget reached; width %d needs %llu extra bytes, "
                    "%llu allocated, budget %llu; current position incomplete\n", n,
                    (unsigned long long)state_counts[n], (unsigned long long)allocated_bytes,
                    (unsigned long long)byte_budget);
            exit(3);
        }
        dense[n] = calloc((size_t)state_counts[n], 1);
        if (!dense[n]) { fprintf(stderr, "strip3_language: allocation failed\n"); exit(3); }
        allocated_bytes += state_counts[n];
    }
    if (dense[n][rank]) return dense[n][rank] - 1;
    int nodes = order(m);
    if (language_flags & 1) {
        int maxdegree = 0, endpoints = 0;
        for (Mask x = m; x; x &= x - 1) {
            int degree = order(neighbours((Mask)1 << lowbit(x)) & m);
            if (degree > maxdegree) maxdegree = degree;
            if (degree < 2) endpoints++;
            if (maxdegree > 2) break;
        }
        if (maxdegree <= 2) {
            int g = (endpoints || nodes == 1) ? path_g[nodes] : (path_g[nodes-3] == 0);
            if (endpoints || nodes == 1) tst.path_hits++; else tst.cycle_hits++;
            save_language(n, rank, g);
            return g;
        }
    }
    Mask vertical = flip_rows(m), both = reverse_bits(m, 3*n), horizontal = flip_rows(both);
    if ((language_flags & 2) && m == both) {
        int safe = 1;
        for (Mask x = m; x && safe; x &= x - 1) {
            int v = lowbit(x), partner = 3*n - 1 - v;
            if (partner == v || (neighbours((Mask)1 << v) & ((Mask)1 << partner))) safe = 0;
        }
        if (safe) { tst.rotation_hits++; save_language(n, rank, 0); return 0; }
    }
    uint64_t seen[2] = {0, 0};
    for (Mask x = m; x; x &= x - 1) {
        int v = lowbit(x), c = v / 3, r = v % 3;
        if (m == vertical && 3*c+2-r < v) continue;
        if (m == horizontal && 3*(n-1-c)+r < v) continue;
        if (m == both && 3*n-1-v < v) continue;
        tst.moves++;
        Mask bit = (Mask)1 << v;
        int g = language_value(m & ~(bit | neighbours(bit)));
        if (g < 128) seen[g >> 6] |= 1ULL << (g & 63);
    }
    int g = 0;
    while (g < 128 && ((seen[g >> 6] >> (g & 63)) & 1)) g++;
    if (g >= 128) { fprintf(stderr, "strip3_language: mex overflow\n"); exit(1); }
    save_language(n, rank, g);
    return g;
}

static int language_value(Mask m) {
    Mask parts[MAXCELLS];
    int count = 0;
    while (m) {
        Mask comp = (Mask)1 << lowbit(m), frontier = comp;
        while (frontier) {
            frontier = neighbours(frontier) & m & ~comp;
            comp |= frontier;
        }
        m &= ~comp;
        Mask key = language_canonical(comp);
        int duplicate = -1;
        if (language_flags & 4)
            for (int i = 0; i < count; i++) if (parts[i] == key) { duplicate = i; break; }
        if (duplicate >= 0) { parts[duplicate] = parts[--count]; tst.cancelled++; }
        else parts[count++] = key;
    }
    int g = 0;
    for (int i = 0; i < count; i++) g ^= language_component(parts[i]);
    return g;
}

static uint64_t read_setting(const char *name, uint64_t fallback, uint64_t maximum) {
    const char *env = getenv(name);
    if (!env) return fallback;
    char *end;
    errno = 0;
    unsigned long long number = strtoull(env, &end, 10);
    if (!*env || *env == '-' || *end || errno || number > maximum) {
        fprintf(stderr, "strip3_language: invalid %s\n", name); exit(2);
    }
    return number;
}

static void language_init(void) {
    for (int c = 0; c < MAXCOLS; c++)
        for (int r = 0; r < 3; r++) rowmask[r] |= (Mask)1 << (3*c+r);
    for (int size = 1; size <= MAXCELLS; size++) {
        uint64_t seen[2] = {0, 0};
        for (int i = 0; i < size; i++) {
            int left = i > 0 ? i - 1 : 0, right = size - i > 1 ? size - i - 2 : 0;
            int g = path_g[left] ^ path_g[right];
            if (g < 128) seen[g >> 6] |= 1ULL << (g & 63);
        }
        while ((seen[path_g[size] >> 6] >> (path_g[size] & 63)) & 1) path_g[size]++;
    }
    state_counts[1] = 6; state_counts[2] = 24;
    for (int n = 3; n <= MAXCOLS; n++) state_counts[n] = completion(n, 0, 0, 4, 0);
    byte_budget = read_setting("STRIP3_LANGUAGE_BYTES", byte_budget, SIZE_MAX);
    language_flags = (int)read_setting("STRIP3_LANGUAGE_REDUCTIONS", 7, 7);
    started = seconds();
}

int main(int argc, char **argv) {
    if (argc < 2) goto usage;
    language_init();
    if ((!strcmp(argv[1], "masks") || !strcmp(argv[1], "ranks")) && argc == 2) {
        char line[128], hex[40], extra;
        int n;
        while (fgets(line, sizeof line, stdin)) {
            if (sscanf(line, "%d %39s %c", &n, hex, &extra) != 2 || n < 1 || n > MAXCOLS) goto usage;
            if (!*hex || strlen(hex) > 32) goto usage;
            Mask m = 0;
            for (char *p = hex; *p; p++) {
                int digit = *p >= '0' && *p <= '9' ? *p-'0' :
                            (*p >= 'a' && *p <= 'f' ? *p-'a'+10 :
                            (*p >= 'A' && *p <= 'F' ? *p-'A'+10 : -1));
                if (digit < 0) goto usage;
                m = (m << 4) | digit;
            }
            if (m >> (3*n)) goto usage;
            if (!strcmp(argv[1], "ranks")) {
                if (!m || !(m & 7)) invalid_mask(m);
                Mask reached = (Mask)1 << lowbit(m), frontier = reached;
                while (frontier) { frontier = neighbours(frontier) & m & ~reached; reached |= frontier; }
                if (reached != m) invalid_mask(m);
                printf("%llu\n", (unsigned long long)language_rank(m));
            } else printf("%d\n", language_value(m));
            fflush(stdout);
        }
    } else if (argc == 3 && (!strcmp(argv[1], "board") || !strcmp(argv[1], "strip"))) {
        char *end;
        long n = strtol(argv[2], &end, 10);
        if (!*argv[2] || *end || n < 1 || n > MAXCOLS) goto usage;
        for (int w = !strcmp(argv[1], "strip") ? 1 : (int)n; w <= n; w++) {
            int g = language_value(full_board(w));
            printf("%d\t%d\n", w, g); fflush(stdout);
        }
    } else goto usage;
    fprintf(stderr, "language: time=%.3f memo=%llu sg_bytes=%llu calls=%llu moves=%llu\n",
            seconds()-started, (unsigned long long)language_entries,
            (unsigned long long)allocated_bytes, (unsigned long long)tst.calls,
            (unsigned long long)tst.moves);
    for (int n = 1; n <= MAXCOLS; n++) free(dense[n]);
    return 0;
usage:
    fprintf(stderr, "usage: strip3_language board N | strip N | masks | ranks\n");
    return 2;
}
