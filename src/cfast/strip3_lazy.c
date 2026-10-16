/* Resumable target queries with bounded work on each unresolved child.
 * A suspended query returns UNKNOWN, never false. Already disproved graph
 * options are retained per (component,target), and finished facts use the
 * ordinary shared component memo. No tuple proof table or stored edges.
 * After an incomplete pass, child quanta grow geometrically; a cursor rotates
 * across options so a hard early branch does not monopolise every retry.
 * Exact SG, N<=42; this is a screening prototype, not a solver for 3x100.
 */
#define TARGET_ENTRY_POINT target_cli_main
#include "strip3_target.c"

enum { UNKNOWN = -1 };
struct pending_query {
    Mask key, done;
    unsigned char k, level, cursor;
};
static struct pending_query *pending;
static size_t pending_capacity, pending_count;
static uint64_t lazy_work, lazy_unknown, lazy_rounds, lazy_skipped, lazy_cached;
static uint64_t work_limit;
static unsigned initial_quantum = 16, growth_bits = 2;
static int lazy_reported;

static void lazy_stats(void) {
    if (lazy_reported) return;
    fprintf(stderr, "lazy work=%llu unknown=%llu rounds=%llu skipped=%llu cached=%llu pending=%zu pending_bytes=%zu quantum=%u growth_bits=%u\n",
            (unsigned long long)lazy_work, (unsigned long long)lazy_unknown,
            (unsigned long long)lazy_rounds, (unsigned long long)lazy_skipped,
            (unsigned long long)lazy_cached, pending_count,
            pending_capacity*sizeof *pending, initial_quantum, growth_bits);
    lazy_reported = 1;
}

static void incomplete(const char *reason) {
    fprintf(stderr, "lazy: %s; incomplete; memo=%zu queries=%llu moves=%llu option_probes=%llu\n", reason, query_entries,
            (unsigned long long)queries, (unsigned long long)query_moves, (unsigned long long)option_probes);
    exit(3);
}

static struct pending_query *pending_slot(Mask m, int k, int create) {
    size_t i = (hash(m)^((uint64_t)k*0x9e3779b97f4a7c15ULL))&(pending_capacity-1);
    while (pending[i].key && (pending[i].key != m || pending[i].k != k)) i = (i+1)&(pending_capacity-1);
    if (!pending[i].key) {
        if (!create) return NULL;
        if (pending_count >= pending_capacity*7/10) incomplete("pending-query budget reached");
        pending[i].key = m; pending[i].k = (unsigned char)k; pending_count++;
    }
    return &pending[i];
}

static uint64_t child_end(uint64_t end, uint64_t quantum) {
    return end-lazy_work < quantum ? end : lazy_work+quantum;
}

static int spend(uint64_t end) {
    if (work_limit && lazy_work >= work_limit) incomplete("work budget reached");
    if (lazy_work >= end) return 0;
    lazy_work++; return 1;
}

static int cached_component(Mask m, int k, int *g) {
    struct query_entry *entry = query_slot(m);
    *g = entry->exact;
    if (*g >= 0) return *g == k;
    if (k > order(m) || (entry->excluded >> k & 1)) return 0;
    return UNKNOWN;
}

static int quick_component(Mask m, struct query_entry *entry) {
    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask both = reverse_bits(m, 3*n) >> (3-h);
    if (m == both) {
        int safe = 1;
        for (Mask x = m; x && safe; x &= x-1) {
            int v = lowbit(x), partner = 3*(n-1-v/3)+h-1-v%3;
            if (partner == v || (neighbours((Mask)1 << v) & ((Mask)1 << partner))) safe = 0;
        }
        if (safe) { set_exact(entry, 0); return 0; }
    }
    int maximum = 0, endpoints = 0, nodes = order(m);
    for (Mask x = m; x; x &= x-1) {
        int degree = order(neighbours((Mask)1 << lowbit(x)) & m);
        if (degree > maximum) maximum = degree;
        if (degree < 2) endpoints++;
        if (maximum > 2) break;
    }
    if (maximum <= 2) {
        int g = endpoints || nodes == 1 ? path_g[nodes] : (path_g[nodes-3] == 0);
        set_exact(entry, g); return g;
    }
    return -1;
}

static int attempt_component(Mask m, int k, uint64_t end);

static int attempt_exact(Mask m, uint64_t end) {
    struct query_entry *entry = query_slot(m);
    if (entry->exact >= 0) return entry->exact;
    for (int k = 0; k <= order(m); k++) {
        int result = attempt_component(m, k, end);
        if (result == UNKNOWN) return UNKNOWN;
        if (result) return k;
        if (entry->exact >= 0) return entry->exact;
    }
    incomplete("inconsistent mex"); return UNKNOWN;
}

/* Return equality and, separately, the exact xor if it became known.
 * No UNKNOWN target may enter the ordinary exact_value loop. */
static int attempt_sum(Mask m, int k, uint64_t end, int *value_out) {
    *value_out = -1;
    if (!spend(end)) return UNKNOWN;
    Mask parts[MAXCELLS]; int count = split_position(m, parts), offset = 0;
    if (!count) { *value_out = 0; return k == 0; }
    for (int i = 1; i < count; i++) {
        int g = attempt_exact(parts[i], end);
        if (g == UNKNOWN) return UNKNOWN;
        offset ^= g;
    }
    int result = attempt_component(parts[0], k^offset, end);
    struct query_entry *entry = query_slot(parts[0]);
    if (entry->exact >= 0) *value_out = entry->exact^offset;
    return result;
}

static int attempt_impl(Mask m, int k, uint64_t end) {
    queries++;
    int g, result = cached_component(m, k, &g);
    if (result != UNKNOWN) { lazy_cached++; return result; }
    if (!spend(end)) return UNKNOWN;
    struct query_entry *entry = query_slot(m);
    g = quick_component(m, entry);
    if (g >= 0) return g == k;
    if (cache_first) {
        int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
        Mask vertical = flip_rows(m) >> (3-h), both = reverse_bits(m, 3*n) >> (3-h);
        Mask horizontal = flip_rows(both) >> (3-h);
        for (int step = 0; step < n; step++) {
            int c = query_order >= 2 ? n/2+(step%2 ? -(step+1)/2 : step/2) : step;
            for (int r = 0; r < h; r++) {
                int v = 3*c+r; Mask bit = (Mask)1 << v;
                if (!(m & bit)) continue;
                if (m == vertical && 3*c+h-1-r < v) continue;
                if (m == horizontal && 3*(n-1-c)+r < v) continue;
                if (m == both && 3*(n-1-c)+h-1-r < v) continue;
                if (!spend(end)) return UNKNOWN;
                option_probes++;
                int child_g = cached_sum(m & ~(bit | neighbours(bit)));
                if (child_g >= 0) {
                    entry->excluded |= (Mask)1 << child_g;
                    if (child_g == k) return 0;
                }
            }
        }
    }
    struct pending_query local = {.key=m, .k=(unsigned char)k};
    struct pending_query *p = pending_slot(m, k, 0);
    if (!p) p = &local;
    unsigned shift = p->level*growth_bits;
    uint64_t quantum = shift >= 62 || initial_quantum > ((UINT64_MAX/4) >> shift)
                           ? UINT64_MAX/4 : (uint64_t)initial_quantum << shift;
    lazy_rounds++;
    int unknown = 0;
    if (k >= frontier_threshold) for (int j = 0; j < k; j++) {
        result = attempt_component(m, j, child_end(end, quantum));
        if (result == 1) return 0;
        if (entry->exact >= 0) return entry->exact == k;
        unknown |= result == UNKNOWN;
        if (lazy_work >= end) goto deferred;
    }
    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask vertical = flip_rows(m) >> (3-h), both = reverse_bits(m, 3*n) >> (3-h);
    Mask horizontal = flip_rows(both) >> (3-h);
    int moves[MAXCELLS], count = 0;
    for (int pass = 0; pass < (query_order ? 2 : 1); pass++) for (int step = 0; step < n; step++) {
        int c = query_order >= 2 ? n/2+(step%2 ? -(step+1)/2 : step/2) : step;
        for (int r = 0; r < h; r++) {
            int v = 3*c+r;
            if (!(m >> v & 1)) continue;
            int prefer_middle = k < 4 || (query_order != 3 && !(query_order == 4 && order(m) == 3*n));
            if (query_order && ((r == 1) != (pass == (prefer_middle ? 0 : 1)))) continue;
            if (m == vertical && 3*c+h-1-r < v) continue;
            if (m == horizontal && 3*(n-1-c)+r < v) continue;
            if (m == both && 3*(n-1-c)+h-1-r < v) continue;
            moves[count++] = v;
        }
    }
    int start = p->cursor%count;
    for (int i = 0; i < count; i++) {
        int index = (start+i)%count, v = moves[index]; Mask bit = (Mask)1 << v;
        if (p->done & bit) { lazy_skipped++; continue; }
        p->cursor = (unsigned char)((index+1)%count);
        if (lazy_work >= end) goto deferred;
        if (move_limit && query_moves >= move_limit) query_budget();
        query_moves++; active_query->move = v;
        Mask child = m & ~(bit | neighbours(bit));
        result = attempt_sum(child, k, child_end(end, quantum), &g);
        if (g >= 0) entry->excluded |= (Mask)1 << g;
        if (result == 1) {
            entry->excluded |= (Mask)1 << k; refutations++; return 0;
        }
        if (result == 0) p->done |= bit;
        else unknown = 1;
    }
    for (int j = 0; j < k; j++) {
        result = attempt_component(m, j, child_end(end, quantum));
        if (result == 1) return 0;
        if (entry->exact >= 0) return entry->exact == k;
        unknown |= result == UNKNOWN;
        if (lazy_work >= end) goto deferred;
    }
    if (!unknown) { set_exact(entry, k); return 1; }
deferred:
    if (p->level < 31) p->level++;
    if (p == &local) *pending_slot(m, k, 1) = local;
    lazy_unknown++;
    return UNKNOWN;
}

static int attempt_component(Mask m, int k, uint64_t end) {
    struct query_frame frame = {m, k, -1, active_query};
    active_query = &frame;
    int result = attempt_impl(m, k, end);
    active_query = frame.parent;
    return result;
}

static int lazy_single(Mask m, int k) {
    int result;
    do { result = attempt_component(m, k, UINT64_MAX); } while (result == UNKNOWN);
    return result;
}

static int configure_lazy(void) {
    const char *env = getenv("STRIP3_LAZY_QUANTUM");
    int quantum = env ? parse_number(env, 1, 1048576) : 16;
    env = getenv("STRIP3_LAZY_GROWTH_BITS");
    int growth = env ? parse_number(env, 1, 4) : 2;
    env = getenv("STRIP3_LAZY_BITS");
    int bits = env ? parse_number(env, 8, 22) : 19;
    env = getenv("STRIP3_LAZY_WORK_LIMIT");
    if (env) {
        char *end;
        if (*env == '-') return 2;
        work_limit = strtoull(env, &end, 10);
        if (!*env || *end || !work_limit) return 2;
    }
    if (quantum < 0 || growth < 0 || bits < 0) return 2;
    initial_quantum = (unsigned)quantum; growth_bits = (unsigned)growth;
    pending_capacity = (size_t)1 << bits;
    pending = calloc(pending_capacity, sizeof *pending);
    if (!pending) return 3;
    atexit(lazy_stats);
    share_options = 1; frontier_threshold = 4;
    component_equality_hook = lazy_single;
    return 0;
}

#ifndef STRIP3_LAZY_EXTERNAL_MAIN
int main(int argc, char **argv) {
    int status = configure_lazy();
    if (status) { fprintf(stderr, "lazy: invalid configuration or allocation failed\n"); return status; }
    status = target_cli_main(argc, argv);
    lazy_stats(); free(pending);
    return status;
}
#endif
