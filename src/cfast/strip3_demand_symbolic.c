/* Demand-driven proof generalisation, not full-domain retrograde iteration.
 * C[k] proves SG=k; W[k] proves SG!=k.  A concrete query selects obligations,
 * but each completed proof learns a Boolean region of induced positions.
 * In the history backend, facts outside the independent-set domain are don't
 * cares and are removed on export. Moves preserve that domain.
 * W is extended only by a legal move to C[k] or a smaller C[j].  C is extended
 * by the universal precondition of the already proved W regions.  Thus no
 * guessed value, unproved state or recursive assumption enters a learned set.
 * Supports 1..100 columns geometrically; scalability is to be measured.
 */
#define STRIP3_REACHABLE_SYMBOLIC
#define STRIP3_PARTIAL_KERNEL_SYMBOLIC
#define main unused_retrograde_main
#include "strip3_symbolic.c"
#undef main

struct position { uint64_t word[5]; };
static struct position closed[300];
static uint64_t demand_queries, demand_moves, demand_hits, demand_p, demand_w;
static size_t gc_watermark;
static int generalise = 2;
static int trace_progress;
struct cube { struct position yes, no; };
static void demand_stats(void);

static int present(struct position m, unsigned v) {
    return (int)((m.word[v/64] >> (v%64)) & 1);
}

static int legal_move(struct position m, unsigned v) {
#ifdef STRIP3_HISTORY_SYMBOLIC
    for (unsigned w = 0; w < 5; w++) if (m.word[w] & closed[v].word[w]) return 0;
    return 1;
#else
    return present(m, v);
#endif
}

static struct position after_move(struct position m, unsigned v) {
#ifdef STRIP3_HISTORY_SYMBOLIC
    m.word[v/64] |= UINT64_C(1) << (v%64);
#else
    for (unsigned w = 0; w < 5; w++) m.word[w] &= ~closed[v].word[w];
#endif
    return m;
}

static Ref legal_formula(unsigned v) {
#ifdef STRIP3_HISTORY_SYMBOLIC
    Ref f = 1;
    for (unsigned u = 0; u < variables; u++) if (present(closed[v], u)) f = bdd_and(f, xvars[u]^1);
    return f;
#else
    return xvars[v];
#endif
}

static int member(Ref f, struct position m) {
    while (f >= 2) {
        struct node node = nodes[f >> 1];
        f = (present(m, node.var) ? node.hi : node.lo) ^ (f & 1);
    }
    return (int)f;
}

#ifdef STRIP3_HISTORY_SYMBOLIC
static Ref set_chosen(Ref f) {
    if (f < 2 || level(f) > current_vertex) return f;
    unsigned id = f >> 1;
    if (restrict_generation[id] == generation) return restrict_cache[id] ^ (f&1);
    struct node node = nodes[id];
    Ref result;
    if (node.var == current_vertex) result = node.hi;
    else {
        Ref lo = set_chosen(node.lo), hi = set_chosen(node.hi);
        result = make_node(node.var, lo, hi);
    }
    restrict_generation[id] = generation; restrict_cache[id] = result;
    return result ^ (f&1);
}
#endif

static Ref inverse_move(Ref f, unsigned v) {
    current_vertex = v;
    if (!++generation) {
        memset(restrict_generation, 0, node_limit*sizeof *restrict_generation);
        generation++;
    }
#ifdef STRIP3_HISTORY_SYMBOLIC
    return set_chosen(f);
#else
    return restrict_zero(f);
#endif
}

static Ref singleton(struct position m) {
    Ref f = 1;
    for (int v = (int)variables-1; v >= 0; v--)
        f = present(m, (unsigned)v) ? make_node((unsigned)v, 0, f) : make_node((unsigned)v, f, 0);
    return f;
}

/* A path to true is an implicant, not an approximation of a false branch.
 * Kill substitutions can discard negative literals at deleted vertices. */
static void explain(Ref f, struct position m, struct cube *cube, int move) {
    if (!member(f, m)) { fprintf(stderr, "demand: false explanation; incomplete\n"); exit(4); }
    while (f >= 2) {
        struct node node = nodes[f >> 1];
        int high = present(m, node.var);
#ifdef STRIP3_HISTORY_SYMBOLIC
        int substituted = move >= 0 && node.var == (unsigned)move;
        if (substituted && !high) { fprintf(stderr, "demand: selected negative literal; incomplete\n"); exit(4); }
#else
        int substituted = move >= 0 && present(closed[move], node.var);
        if (substituted && high) { fprintf(stderr, "demand: killed positive literal; incomplete\n"); exit(4); }
#endif
        if (!substituted) {
            struct position *literals = high ? &cube->yes : &cube->no;
            literals->word[node.var/64] |= UINT64_C(1) << (node.var%64);
        }
        f = (high ? node.hi : node.lo) ^ (f & 1);
    }
}

static Ref cube_ref(struct cube cube) {
    for (unsigned w = 0; w < 5; w++) if (cube.yes.word[w] & cube.no.word[w]) {
        fprintf(stderr, "demand: inconsistent explanation; incomplete\n"); exit(4);
    }
    Ref f = 1;
    for (int v = (int)variables-1; v >= 0; v--) {
        if (present(cube.yes, (unsigned)v)) f = make_node((unsigned)v, 0, f);
        else if (present(cube.no, (unsigned)v)) f = make_node((unsigned)v, f, 0);
    }
    return f;
}

/* Query stacks contain only concrete positions and integers, never BDD refs.
 * Collection here can therefore remap every live root, including all W[k]. */
static void demand_collect(void) {
    if (gc_enabled && node_count >= gc_watermark) {
        Ref a = 0, b = 0, c = 0;
        collect_live(&a, &b, &c);
        gc_watermark = node_count*2;
        if (gc_watermark < node_limit*3/4) gc_watermark = node_limit*3/4;
        if (gc_watermark >= node_limit) gc_watermark = node_count+(node_limit-node_count)/2;
    }
}

static void learn_w(unsigned k, Ref f) {
    partial_winning[k] = bdd_or(partial_winning[k], f);
    demand_w++;
}

static int demand_equal(struct position m, unsigned k) {
    demand_queries++;
    if (trace_progress && demand_queries >= 256 && !(demand_queries & (demand_queries-1))) demand_stats();
    demand_collect();
    if (member(certificate_sets[k], m)) { demand_hits++; return 1; }
    if (member(partial_winning[k], m)) { demand_hits++; return 0; }
    for (unsigned j = 0; j < k; j++) if (member(certificate_sets[j], m)) {
        learn_w(k, certificate_sets[j]);
        demand_hits++; return 0;
    }
    unsigned n = variables/3;
    for (unsigned pass = 0; pass < 2; pass++) for (unsigned step = 0; step < n; step++) {
        int c = (int)n/2 + (step%2 ? -((int)step+1)/2 : (int)step/2);
        for (unsigned r = 0; r < 3; r++) {
            unsigned v = 3*(unsigned)c+r;
            if ((r == 1) != (pass == 0) || !legal_move(m, v)) continue;
            struct position next = after_move(m, v);
            demand_moves++;
            if (demand_equal(next, k)) {
                Ref f;
                if (generalise == 2) {
                    struct cube cube = {{{0}}, {{0}}};
#ifdef STRIP3_HISTORY_SYMBOLIC
                    cube.no = closed[v];
#else
                    cube.yes.word[v/64] |= UINT64_C(1) << (v%64);
#endif
                    explain(certificate_sets[k], next, &cube, (int)v);
                    f = cube_ref(cube);
                } else if (generalise) {
                    Ref inverse = inverse_move(certificate_sets[k], v);
                    f = bdd_and(legal_formula(v), inverse);
                } else f = singleton(m);
                learn_w(k, f);
                return 0;
            }
        }
    }
    for (unsigned j = 0; j < k; j++) if (demand_equal(m, j)) {
        learn_w(k, certificate_sets[j]);
        return 0;
    }
    /* Every legal concrete child is now in W[k].  Generalisation checks ALL
     * vertices, allowing extra legal moves only where W[k] already covers
     * their result.  Smaller-Nim obligations use the proved W[j] regions. */
    Ref region = 1;
    if (generalise == 2) {
        struct cube cube = {{{0}}, {{0}}};
        for (unsigned v = 0; v < variables; v++) {
            if (!legal_move(m, v)) {
#ifdef STRIP3_HISTORY_SYMBOLIC
                for (unsigned u = 0; u < variables; u++) if (present(m, u) && present(closed[v], u)) {
                    cube.yes.word[u/64] |= UINT64_C(1) << (u%64);
                    break;
                }
#else
                cube.no.word[v/64] |= UINT64_C(1) << (v%64);
#endif
            }
            else {
                struct position next = after_move(m, v);
                explain(partial_winning[k], next, &cube, (int)v);
            }
        }
        for (unsigned j = 0; j < k; j++) explain(partial_winning[j], m, &cube, -1);
        region = cube_ref(cube);
    } else if (generalise) {
        for (unsigned v = 0; v < variables; v++) {
            Ref inverse = inverse_move(partial_winning[k], v);
            region = bdd_and(region, bdd_or(legal_formula(v)^1, inverse));
        }
        for (unsigned j = 0; j < k; j++) region = bdd_and(region, partial_winning[j]);
    } else region = singleton(m);
    if (!member(region, m)) { fprintf(stderr, "demand: uncovered completed proof; incomplete\n"); exit(4); }
    certificate_sets[k] = bdd_or(certificate_sets[k], region);
    demand_p++;
    return 1;
}

static int demand_value(struct position m) {
    for (unsigned k = 0; k <= variables; k++) if (demand_equal(m, k)) return (int)k;
    fprintf(stderr, "demand: invalid mex; incomplete\n"); exit(4);
}

static void demand_stats(void) {
    fprintf(stderr, "demand queries=%llu moves=%llu hits=%llu p_facts=%llu w_facts=%llu created_nodes=%llu peak_nodes=%zu retained_nodes=%zu gc=%llu apply_calls=%llu generalise=%d\n",
            (unsigned long long)demand_queries, (unsigned long long)demand_moves,
            (unsigned long long)demand_hits, (unsigned long long)demand_p,
            (unsigned long long)demand_w, (unsigned long long)created_nodes,
            peak_nodes, node_count, (unsigned long long)collections,
            (unsigned long long)apply_calls, generalise);
}

int main(int argc, char **argv) {
    if (argc == 5 && !strcmp(argv[1], "board") && !strcmp(argv[3], "--certificate")) {
        certificate_path = argv[4]; argc = 3;
    }
    int position_query = argc == 4 && !strcmp(argv[1], "position");
    if (!position_query && (argc != 3 || (strcmp(argv[1], "board") && strcmp(argv[1], "table")))) return 2;
    char *end;
    long n = strtol(argv[2], &end, 10);
    int table = !strcmp(argv[1], "table");
    if (!*argv[2] || *end || n < 1 || n > 100 || (table && n > 6)) return 2;
    variables = (unsigned)(3*n);
    const char *env = getenv("STRIP3_SYMBOLIC_BITS");
    long bits = env ? strtol(env, &end, 10) : 20;
    if (bits < 8 || bits > 23 || (env && (!*env || *end))) return 2;
    env = getenv("STRIP3_DEMAND_GENERALISE");
    if (env && strcmp(env, "0") && strcmp(env, "1") && strcmp(env, "2")) return 2;
    generalise = env ? atoi(env) : 2;
    env = getenv("STRIP3_SYMBOLIC_GC");
    if (env && strcmp(env, "0") && strcmp(env, "1")) return 2;
    gc_enabled = !env || !strcmp(env, "1");
    env = getenv("STRIP3_DEMAND_TRACE");
    if (env && strcmp(env, "0") && strcmp(env, "1")) return 2;
    trace_progress = env && !strcmp(env, "1");
    node_limit = (size_t)1 << bits; unique_cap = node_limit*2; cache_cap = node_limit;
    gc_watermark = node_limit*3/4;
    nodes = calloc(node_limit, sizeof *nodes);
    unique = calloc(unique_cap, sizeof *unique);
    cache = calloc(cache_cap, sizeof *cache);
    restrict_generation = calloc(node_limit, sizeof *restrict_generation);
    restrict_cache = calloc(node_limit, sizeof *restrict_cache);
    compact_map = calloc(node_limit, sizeof *compact_map);
    marked = calloc(node_limit, sizeof *marked);
    if (!nodes || !unique || !cache || !restrict_generation || !restrict_cache || !compact_map || !marked) return 3;
    atexit(demand_stats);
    certificate_count = variables+1;
    for (unsigned v = 0; v < variables; v++) {
        xvars[v] = make_node(v, 0, 1);
        for (unsigned u = 0; u < variables; u++) {
            int dc = (int)(v/3)-(int)(u/3), dr = (int)(v%3)-(int)(u%3);
            if (abs(dc)+abs(dr) <= 1) closed[v].word[u/64] |= UINT64_C(1) << (u%64);
        }
    }
    struct position full = {{0}}, empty = {{0}};
#ifdef STRIP3_HISTORY_SYMBOLIC
    (void)empty;
    domain = 1;
    /* Prepend local clauses to the existing suffix. Appending them left to
     * right retains quadratically many obsolete prefix diagrams at setup. */
    for (int v = (int)variables-1; v >= 0; v--) for (unsigned u = (unsigned)v+1; u < variables; u++)
        if (present(closed[v], u)) domain = bdd_and(domain, bdd_or(xvars[v]^1, xvars[u]^1));
    Ref terminal = 1;
    for (int v = (int)variables-1; v >= 0; v--) terminal = bdd_and(terminal, legal_formula((unsigned)v)^1);
    for (unsigned k = 0; k < certificate_count; k++) certificate_sets[k] = partial_winning[k] = domain^1;
    certificate_sets[0] = bdd_or(certificate_sets[0], terminal);
#else
    for (unsigned v = 0; v < variables; v++) full.word[v/64] |= UINT64_C(1) << (v%64);
    certificate_sets[0] = singleton(empty);
#endif
    if (position_query) {
        memset(&full, 0, sizeof full);
        size_t length = strlen(argv[3]);
        if (!length || length > 75) return 2;
        for (size_t i = 0; i < length; i++) {
            char c = argv[3][length-1-i];
            unsigned digit;
            if (c >= '0' && c <= '9') digit = (unsigned)(c-'0');
            else if (c >= 'a' && c <= 'f') digit = (unsigned)(c-'a')+10;
            else if (c >= 'A' && c <= 'F') digit = (unsigned)(c-'A')+10;
            else return 2;
            for (unsigned bit = 0; bit < 4; bit++) if ((digit >> bit) & 1) {
                size_t v = 4*i+bit;
                if (v >= variables) return 2;
                full.word[v/64] |= UINT64_C(1) << (v%64);
            }
        }
#ifdef STRIP3_HISTORY_SYMBOLIC
        if (!member(domain, full)) return 2;
#endif
    }
    int g;
    if (table) {
        uint64_t count = UINT64_C(1) << variables;
        int *values = malloc((size_t)count*sizeof *values);
        if (!values) return 3;
        for (uint64_t mask = 0; mask < count; mask++) {
            struct position m = {{mask, 0, 0, 0, 0}};
#ifdef STRIP3_HISTORY_SYMBOLIC
            if (!member(domain, m)) { values[mask] = -1; continue; }
#endif
            values[mask] = demand_value(m);
        }
        for (uint64_t mask = 0; mask < count; mask++) printf("%d\n", values[mask]);
        free(values);
    } else {
        g = demand_value(full);
        if (certificate_path) {
            /* Restriction to a move-closed R preserves the partial-kernel
             * conditions.  Domain construction is paid only for export. */
            demand_collect();
#ifndef STRIP3_HISTORY_SYMBOLIC
            domain = reachable_domain((unsigned)n);
#endif
            for (int k = 0; k <= g; k++) certificate_sets[k] = bdd_and(domain, certificate_sets[k]);
            if (!write_kernel_certificate((unsigned)n, (unsigned)g, certificate_sets[g])) {
                fprintf(stderr, "demand: certificate write failed; incomplete\n"); return 3;
            }
        }
        if (position_query) printf("position\t%ld\t%s\t%d\n", n, argv[3], g);
        else printf("3\t%ld\t%d\n", n, g);
    }
    return 0;
}
