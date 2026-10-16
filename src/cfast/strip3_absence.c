/* Mex bounds from absent option values, before determining the exact mex.
 * A(G,j) = no legal option of G has SG j implies g(G) <= j; it does NOT
 * imply g(G)=j.  Reuse this upper bound to refute every larger target.
 * Actual option witnesses (the inherited excluded bitset) give a lower bound.
 * Equality requires both bounds to meet.  All recursive equality calls in this
 * implementation delete vertices, including the calls used to fill lower gaps.
 * Screening prototype: same 42-column limitation as strip3_target.c.
 */
#define TARGET_ENTRY_POINT target_cli_main
#include "strip3_target.c"

static unsigned char *upper_codes; /* 0 = |V|; otherwise upper bound + 1 */
static int absence_probe = 1;     /* 0=none, 1=k-1, 2=k/2 */
static uint64_t absence_calls, absence_cached, bound_refutations, bound_updates;
static uint64_t bound_cache_refutations;
static int stats_reported;

static void report_absence(void) {
    if (stats_reported) return;
    fprintf(stderr, "absence calls=%llu cached=%llu bound_refutations=%llu bound_cache_refutations=%llu bound_updates=%llu probe=%d extra_bytes=%zu\n",
            (unsigned long long)absence_calls, (unsigned long long)absence_cached,
            (unsigned long long)bound_refutations, (unsigned long long)bound_cache_refutations,
            (unsigned long long)bound_updates, absence_probe, upper_codes ? query_capacity : 0);
    stats_reported = 1;
}

static int upper_bound(struct query_entry *e, int nodes) {
    unsigned char code = upper_codes[e-query_table];
    return code ? code-1 : nodes;
}

static int lower_bound(struct query_entry *e) {
    int k = 0;
    while ((e->excluded >> k) & 1) k++;
    return k;
}

static void infer_exact(struct query_entry *e, int nodes) {
    int lower = lower_bound(e), upper = upper_bound(e, nodes);
    if (lower > upper) {
        fprintf(stderr, "absence: inconsistent proven bounds; incomplete\n");
        exit(4);
    }
    if (lower == upper) set_exact(e, lower);
}

static void tighten_upper(struct query_entry *e, int nodes, int k) {
    if (k < upper_bound(e, nodes)) {
        upper_codes[e-query_table] = (unsigned char)(k+1);
        bound_updates++;
    }
    infer_exact(e, nodes);
}

/* Called only on a canonical connected position.  Every bit published in
 * e->excluded is an ACTUAL option value, not an arbitrary negative equality.
 * A cached upper bound is accompanied by absence of that very option value;
 * larger values need not be absent.  |V| is absent because every child loses
 * at least one vertex and has at most |V|-1 legal moves.
 */
static int no_option(Mask m, int k, struct query_entry *e) {
    absence_calls++;
    int nodes = order(m);
    if ((e->excluded >> k) & 1) { absence_cached++; return 0; }
    if (k >= nodes || k == upper_bound(e, nodes)) { absence_cached++; return 1; }
    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask vertical = flip_rows(m) >> (3-h);
    Mask both = reverse_bits(m, 3*n) >> (3-h);
    Mask horizontal = flip_rows(both) >> (3-h);
    for (int pass = 0; pass < (query_order ? 2 : 1); pass++) {
        for (int step = 0; step < n; step++) {
            int c = query_order >= 2 ? n/2 + (step%2 ? -(step+1)/2 : step/2) : step;
            for (int r = 0; r < h; r++) {
                int v = 3*c+r;
                if (!(m >> v & 1)) continue;
                int prefer_middle = k < 4 || (query_order != 3 && !(query_order == 4 && nodes == 3*n));
                if (query_order && ((r == 1) != (pass == (prefer_middle ? 0 : 1)))) continue;
                if (m == vertical && 3*c+h-1-r < v) continue;
                if (m == horizontal && 3*(n-1-c)+r < v) continue;
                if (m == both && 3*(n-1-c)+h-1-r < v) continue;
                if (move_limit && query_moves >= move_limit) query_budget();
                query_moves++;
                active_query->move = v;
                Mask bit = (Mask)1 << v;
                int matched = equality(m & ~(bit | neighbours(bit)), k);
                if (share_options && last_equality_value >= 0)
                    e->excluded |= (Mask)1 << last_equality_value;
                if (matched) {
                    e->excluded |= (Mask)1 << k;
                    refutations++;
                    infer_exact(e, nodes);
                    return 0;
                }
            }
        }
    }
    tighten_upper(e, nodes, k);
    return 1;
}

static int bounded_equality(Mask m, int k) {
    queries++;
    int nodes = order(m);
    if (k > nodes) return 0;
    if (!upper_codes) {
        upper_codes = calloc(query_capacity, sizeof *upper_codes);
        if (!upper_codes) { fprintf(stderr, "absence: allocation failed; incomplete\n"); exit(3); }
    }
    struct query_entry *e = query_slot(m);
    if (e->exact >= 0) return e->exact == k;
    if ((e->excluded >> k) & 1) return 0;
    if (k > upper_bound(e, nodes)) { bound_refutations++; bound_cache_refutations++; return 0; }

    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask both = reverse_bits(m, 3*n) >> (3-h);
    if (m == both) {
        int safe = 1;
        for (Mask x = m; x && safe; x &= x-1) {
            int v = lowbit(x), partner = 3*(n-1-v/3)+h-1-v%3;
            if (partner == v || (neighbours((Mask)1 << v) & ((Mask)1 << partner))) safe = 0;
        }
        if (safe) { set_exact(e, 0); return k == 0; }
    }
    int maxdegree = 0, endpoints = 0;
    for (Mask x = m; x; x &= x-1) {
        int degree = order(neighbours((Mask)1 << lowbit(x)) & m);
        if (degree > maxdegree) maxdegree = degree;
        if (degree < 2) endpoints++;
        if (maxdegree > 2) break;
    }
    if (maxdegree <= 2) {
        int g = endpoints || nodes == 1 ? path_g[nodes] : (path_g[nodes-3] == 0);
        set_exact(e, g); return g == k;
    }
    if (cache_first) {
        Mask vertical = flip_rows(m) >> (3-h);
        Mask horizontal = flip_rows(both) >> (3-h);
        for (int step = 0; step < n; step++) {
            int c = query_order >= 2 ? n/2 + (step%2 ? -(step+1)/2 : step/2) : step;
            for (int r = 0; r < h; r++) {
                int v = 3*c+r;
                if (!(m >> v & 1)) continue;
                if (m == vertical && 3*c+h-1-r < v) continue;
                if (m == horizontal && 3*(n-1-c)+r < v) continue;
                if (m == both && 3*(n-1-c)+h-1-r < v) continue;
                Mask bit = (Mask)1 << v;
                option_probes++;
                int g = cached_sum(m & ~(bit | neighbours(bit)));
                if (g >= 0) {
                    e->excluded |= (Mask)1 << g;
                    if (g == k) { infer_exact(e, nodes); return 0; }
                }
            }
        }
    }
    infer_exact(e, nodes);
    if (e->exact >= 0) return e->exact == k;
    if ((e->excluded >> k) & 1) return 0;
    if (absence_probe && k >= frontier_threshold) {
        int trial = absence_probe == 1 ? k-1 : k/2;
        if (trial >= lower_bound(e) && no_option(m, trial, e)) {
            bound_refutations++; return 0;
        }
        if (e->exact >= 0) return e->exact == k;
    }
    if (!no_option(m, k, e)) return 0;
    if (e->exact >= 0) return e->exact == k;
    if (upper_bound(e, nodes) < k) { bound_refutations++; return 0; }
    for (int j = 0; j < k; j++) if (!(e->excluded >> j & 1)) {
        if (no_option(m, j, e)) return 0;
        if (e->exact >= 0) return e->exact == k;
    }
    infer_exact(e, nodes);
    if (e->exact != k) { fprintf(stderr, "absence: incomplete mex\n"); exit(4); }
    return 1;
}

int main(int argc, char **argv) {
    const char *env = getenv("STRIP3_ABSENCE_PROBE");
    if (env) absence_probe = parse_number(env, 0, 2);
    if (absence_probe < 0) { fprintf(stderr, "absence: invalid probe mode\n"); return 2; }
    share_options = 1;
    frontier_threshold = 4;
    component_equality_hook = bounded_equality;
    atexit(report_absence);
    int status = target_cli_main(argc, argv);
    report_absence();
    free(upper_codes);
    return status;
}
