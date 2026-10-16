/* Experimental target-nimber search, NOT a solver for 3x100.
 * Reuse bit geometry/normalisation, but neither the SG recursion nor its memo.
 * E(G,k) <=> no option has SG k AND no smaller nimber equals SG(G).
 * Test options first: one matching option refutes E without a complete mex.
 * For disconnected positions, evaluate all but the largest component exactly
 * and xor their values into the query for the largest component.
 * CLI: board N | equal N K | masks (N HEX) | queries (N HEX K).
 * N <= 42 is an explicit limitation of this screening prototype.
 */
#define main original_strip3_main
#include "strip3.c"
#undef main

struct query_entry { Mask key, excluded; int exact; };
static struct query_entry *query_table;
static size_t query_capacity, query_entries;
static uint64_t queries, query_moves, exact_entries, refutations;
static int query_order = 2;  /* 0=natural, 1=middle first, 2=middle/centre first */
static int low_first;

static struct query_entry *query_slot(Mask m) {
    size_t i = hash(m) & (query_capacity - 1);
    while (query_table[i].key && query_table[i].key != m)
        i = (i + 1) & (query_capacity - 1);
    if (!query_table[i].key) {
        if (query_entries >= query_capacity * 7 / 10) {
            fprintf(stderr, "target: memo budget reached; query incomplete\n");
            exit(3);
        }
        query_table[i].key = m;
        query_table[i].exact = -1;
        query_entries++;
    }
    return &query_table[i];
}

static int equality(Mask m, int k);
static int exact_value(Mask m);

static int split_position(Mask m, Mask *parts) {
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
        for (int i = 0; i < count; i++) if (parts[i] == key) { duplicate = i; break; }
        if (duplicate >= 0) parts[duplicate] = parts[--count];
        else parts[count++] = key;
    }
    /* Keep the largest component as a query; solve smaller ones exactly. */
    for (int i = 1; i < count; i++) if (order(parts[i]) > order(parts[0])) {
        Mask swap = parts[0]; parts[0] = parts[i]; parts[i] = swap;
    }
    return count;
}

static void set_exact(struct query_entry *e, int g) {
    if (e->exact < 0) exact_entries++;
    e->exact = g;
}

static int equality_component(Mask m, int k) {
    queries++;
    int nodes = order(m);
    if (k > nodes) return 0;  /* mex of at most |V| options */
    struct query_entry *e = query_slot(m);
    if (e->exact >= 0) return e->exact == k;
    if ((e->excluded >> k) & 1) return 0;

    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask vertical = flip_rows(m) >> (3 - h);
    Mask both = reverse_bits(m, 3*n) >> (3 - h);
    Mask horizontal = flip_rows(both) >> (3 - h);
    if (m == both) {
        int safe = 1;
        for (Mask x = m; x && safe; x &= x - 1) {
            int v = lowbit(x), partner = 3*(n-1-v/3) + h-1-v%3;
            if (partner == v || (neighbours((Mask)1 << v) & ((Mask)1 << partner))) safe = 0;
        }
        if (safe) { set_exact(e, 0); return k == 0; }
    }
    int degree_max = 0, endpoints = 0;
    for (Mask x = m; x; x &= x - 1) {
        int degree = order(neighbours((Mask)1 << lowbit(x)) & m);
        if (degree > degree_max) degree_max = degree;
        if (degree < 2) endpoints++;
        if (degree_max > 2) break;
    }
    if (degree_max <= 2) {
        int g = endpoints || nodes == 1 ? path_g[nodes] : (path_g[nodes-3] == 0);
        set_exact(e, g);
        return g == k;
    }

    /* Optional Nim move ordering: try the zero target before graph moves.
     * This can refute a high target without proving absence of that target
     * among all graph options. It is an ablation, not an assumed speedup. */
    if (low_first && k && !(e->excluded & 1) && equality_component(m, 0)) return 0;

    int move_list[MAXCELLS], move_count = 0;
    if (query_order == 2) {
        for (int i = 0; i < n; i++) {
            int c = n/2 + (i % 2 ? -(i+1)/2 : i/2);
            for (int r = 0; r < h; r++) if ((m >> (3*c+r)) & 1)
                move_list[move_count++] = 3*c+r;
        }
    } else for (Mask x = m; x; x &= x - 1) move_list[move_count++] = lowbit(x);
    for (int pass = 0; pass < (query_order ? 2 : 1); pass++) {
        for (int i = 0; i < move_count; i++) {
            int v = move_list[i], c = v / 3, r = v % 3;
            if (query_order && ((r == 1) != (pass == 0))) continue;
            if (m == vertical && 3*c+h-1-r < v) continue;
            if (m == horizontal && 3*(n-1-c)+r < v) continue;
            if (m == both && 3*(n-1-c)+h-1-r < v) continue;
            query_moves++;
            Mask bit = (Mask)1 << v;
            if (equality(m & ~(bit | neighbours(bit)), k)) {
                e->excluded |= (Mask)1 << k;
                refutations++;
                return 0;
            }
        }
    }
    /* No option equals k. It is the mex iff every smaller value is excluded.
     * Same-position recursion decreases k; move recursion decreases |V|. */
    for (int j = 0; j < k; j++) {
        if (equality_component(m, j)) return 0;
        if (e->exact >= 0) return e->exact == k;
    }
    set_exact(e, k);
    return 1;
}

static int equality(Mask m, int k) {
    Mask parts[MAXCELLS];
    int count = split_position(m, parts);
    if (!count) return k == 0;
    for (int i = 1; i < count; i++) k ^= exact_value(parts[i]);
    return equality_component(parts[0], k);
}

static int exact_value(Mask m) {
    Mask parts[MAXCELLS];
    int count = split_position(m, parts), g = 0;
    for (int i = 0; i < count; i++) {
        int k = 0;
        while (!equality_component(parts[i], k)) k++;
        g ^= k;
    }
    return g;
}

static int parse_number(const char *s, int lo, int hi) {
    char *end;
    long n = strtol(s, &end, 10);
    if (!*s || *end || n < lo || n > hi) return -1;
    return (int)n;
}

int main(int argc, char **argv) {
    if (argc < 2) goto usage;
    int board = !strcmp(argv[1], "board"), eq = !strcmp(argv[1], "equal");
    int batch = !strcmp(argv[1], "masks"), qs = !strcmp(argv[1], "queries");
    int n = 0, k = 0;
    if (board || eq) {
        if (argc != (eq ? 4 : 3)) goto usage;
        n = parse_number(argv[2], 1, MAXCOLS);
        if (eq) k = parse_number(argv[3], 0, 127);
        if (n < 0 || k < 0) goto usage;
    } else if ((!batch && !qs) || argc != 2) goto usage;
    const char *env = getenv("STRIP3_TARGET_BITS");
    int bits = env ? parse_number(env, 8, 26) : 20;
    env = getenv("STRIP3_TARGET_ORDER");
    if (env) query_order = parse_number(env, 0, 2);
    env = getenv("STRIP3_TARGET_LOW_FIRST");
    if (env) low_first = parse_number(env, 0, 1);
    if (bits < 0 || query_order < 0 || low_first < 0) goto usage;
    query_capacity = (size_t)1 << bits;
    query_table = calloc(query_capacity, sizeof *query_table);
    if (!query_table) { fprintf(stderr, "target: allocation failed\n"); return 3; }
    for (int c = 0; c < MAXCOLS; c++)
        for (int r = 0; r < 3; r++) rowmask[r] |= (Mask)1 << (3*c+r);
    for (int size = 1; size <= MAXCELLS; size++) {
        Mask seen = 0;
        for (int i = 0; i < size; i++) {
            int left = i > 0 ? i-1 : 0, right = size-i > 1 ? size-i-2 : 0;
            seen |= (Mask)1 << (path_g[left] ^ path_g[right]);
        }
        while ((seen >> path_g[size]) & 1) path_g[size]++;
    }
    double start = seconds();
    if (board) {
        int g = 0, match;
        do {
            match = equality(full_board(n), g);
            fprintf(stderr, "stage k=%d equal=%d time=%.6f memo=%zu queries=%llu moves=%llu\n",
                    g, match, seconds()-start, query_entries,
                    (unsigned long long)queries, (unsigned long long)query_moves);
            fflush(stderr);  /* retain completed obligations even on a later timeout */
            if (!match) g++;
        } while (!match);
        printf("3\t%d\t%d\n", n, g);
    }
    else if (eq) printf("%d\n", equality(full_board(n), k));
    else {
        char line[160], hex[40], extra;
        while (fgets(line, sizeof line, stdin)) {
            int fields = qs ? sscanf(line, "%d %39s %d %c", &n, hex, &k, &extra)
                            : sscanf(line, "%d %39s %c", &n, hex, &extra);
            if (fields != (qs ? 3 : 2) || n < 1 || n > MAXCOLS || strlen(hex) > 32 || k < 0 || k > 127) goto usage;
            Mask m = 0;
            for (const char *p = hex; *p; p++) {
                int digit = *p >= '0' && *p <= '9' ? *p-'0' :
                    *p >= 'a' && *p <= 'f' ? *p-'a'+10 :
                    *p >= 'A' && *p <= 'F' ? *p-'A'+10 : -1;
                if (digit < 0) goto usage;
                m = (m << 4) | digit;
            }
            if (m & ~full_board(n)) goto usage;
            printf("%d\n", qs ? equality(m, k) : exact_value(m));
        }
    }
    fprintf(stderr, "target time=%.6f memo=%zu exact=%llu queries=%llu moves=%llu refutations=%llu bytes_per_slot=%zu order=%d low_first=%d\n",
            seconds()-start, query_entries, (unsigned long long)exact_entries,
            (unsigned long long)queries, (unsigned long long)query_moves,
            (unsigned long long)refutations, sizeof *query_table, query_order, low_first);
    free(query_table);
    return 0;
usage:
    fprintf(stderr, "usage: strip3_target board N | equal N K | masks | queries (N=1..42, K=0..127)\n");
    return 2;
}
