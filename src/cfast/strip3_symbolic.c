/* Experimental symbolic retrograde SG computation using a complemented ROBDD.
 * One Boolean function represents a set of induced masks. It is NOT assumed
 * that the frontier game has finitely many states. All references use exact
 * unique-table keys, and a hard node limit fails without emitting a result.
 * Usage: board N | table N (all 2^(3N) induced masks, small N only).
 */
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint32_t Ref;
struct node { Ref lo, hi; unsigned var; };
struct cache { Ref a, b, value; };
static struct node *nodes;
static uint32_t *unique, *restrict_generation;
static Ref *restrict_cache;
static struct cache *cache;
static size_t node_count = 1, node_limit, unique_cap, cache_cap;
static uint32_t generation;
static unsigned variables, current_vertex;
static Ref xvars[300];
static uint64_t apply_calls;
static Ref domain = 1;
#ifdef STRIP3_REACHABLE_SYMBOLIC
static uint64_t created_nodes, collections;
static size_t peak_nodes;
static uint32_t *compact_map;
static unsigned char *marked;
static int gc_enabled = 1;
static Ref certificate_sets[301];
static unsigned certificate_count;
static const char *certificate_path;
#ifdef STRIP3_PARTIAL_KERNEL_SYMBOLIC
static Ref partial_winning[301];
#endif
#endif

static uint64_t mix(uint64_t h) {
    h ^= h >> 30; h *= 0xbf58476d1ce4e5b9ULL;
    h ^= h >> 27; h *= 0x94d049bb133111ebULL;
    return h ^ (h >> 31);
}

static unsigned level(Ref f) { return f < 2 ? variables : nodes[f >> 1].var; }

static Ref make_node(unsigned var, Ref lo, Ref hi) {
    if (lo == hi) return lo;
    Ref complement = lo & 1;
    lo ^= complement; hi ^= complement;
    size_t index = mix(((uint64_t)lo << 32) ^ hi ^ ((uint64_t)var*0x9e3779b97f4a7c15ULL)) & (unique_cap-1);
    while (unique[index]) {
        uint32_t id = unique[index];
        if (nodes[id].var == var && nodes[id].lo == lo && nodes[id].hi == hi) return (id<<1)^complement;
        index = (index+1) & (unique_cap-1);
    }
    if (node_count == node_limit) {
        fprintf(stderr, "symbolic: node budget reached; incomplete; nodes=%zu apply_calls=%llu\n",
                node_count, (unsigned long long)apply_calls);
#ifdef STRIP3_REACHABLE_SYMBOLIC
        fprintf(stderr, "reachable_stats created_nodes=%llu peak_nodes=%zu retained_nodes=%zu gc=%llu apply_calls=%llu\n",
                (unsigned long long)created_nodes, peak_nodes, node_count, (unsigned long long)collections,
                (unsigned long long)apply_calls);
#endif
        exit(3);
    }
    uint32_t id = (uint32_t)node_count++;
#ifdef STRIP3_REACHABLE_SYMBOLIC
    created_nodes++;
    if (node_count > peak_nodes) peak_nodes = node_count;
#endif
    nodes[id] = (struct node){lo, hi, var};
    unique[index] = id;
    return (id<<1)^complement;
}

static Ref bdd_and(Ref a, Ref b) {
    apply_calls++;
    if (!a || !b || a == (b^1)) return 0;
    if (a == 1 || a == b) return b;
    if (b == 1) return a;
    if (a > b) { Ref temp = a; a = b; b = temp; }
    size_t i = mix(((uint64_t)a << 32) | b) & (cache_cap-1);
    if (cache[i].a == a && cache[i].b == b) return cache[i].value;
    unsigned la = level(a), lb = level(b), var = la < lb ? la : lb;
    Ref al = a, ah = a, bl = b, bh = b;
    if (la == var) { al = nodes[a>>1].lo ^ (a&1); ah = nodes[a>>1].hi ^ (a&1); }
    if (lb == var) { bl = nodes[b>>1].lo ^ (b&1); bh = nodes[b>>1].hi ^ (b&1); }
    Ref low = bdd_and(al, bl), high = bdd_and(ah, bh);
    Ref result = make_node(var, low, high);
    cache[i] = (struct cache){a, b, result};
    return result;
}

static Ref bdd_or(Ref a, Ref b) { return bdd_and(a^1, b^1)^1; }

static int killed(unsigned var) {
    int r = (int)(var%3), c = (int)(var/3);
    int vr = (int)(current_vertex%3), vc = (int)(current_vertex/3);
    return abs(r-vr)+abs(c-vc) <= 1;
}

static Ref restrict_zero(Ref f) {
    if (f < 2) return f;
    uint32_t id = f >> 1;
    if (restrict_generation[id] == generation) return restrict_cache[id] ^ (f&1);
    struct node node = nodes[id];
    Ref value;
    if (killed(node.var)) value = restrict_zero(node.lo);
    else {
        Ref lo = restrict_zero(node.lo), hi = restrict_zero(node.hi);
        value = make_node(node.var, lo, hi);
    }
    restrict_generation[id] = generation;
    restrict_cache[id] = value;
    return value ^ (f&1);
}

/* pre(F)(S) = there exists an available vertex whose move reaches F. */
static Ref predecessor(Ref f) {
    Ref result = 0;
    for (unsigned v = 0; v < variables; v++) {
        current_vertex = v;
        generation++;
        if (!generation) {
            memset(restrict_generation, 0, node_limit*sizeof *restrict_generation);
            generation++;
        }
        Ref legal = domain == 1 ? xvars[v] : bdd_and(domain, xvars[v]);
        result = bdd_or(result, bdd_and(legal, restrict_zero(f)));
    }
    return result;
}

static int evaluate(Ref f, uint64_t mask, int full) {
    while (f >= 2) {
        struct node node = nodes[f>>1];
        int high = full || ((mask >> node.var) & 1);
        f = (high ? node.hi : node.lo) ^ (f&1);
    }
    return f;
}

#ifdef STRIP3_REACHABLE_SYMBOLIC
/* Collection only occurs between fixed-point iterations, when all live roots
 * are explicit. No collection can invalidate a recursive apply/restrict frame. */
static void mark_live(Ref f) {
    if (f < 2 || marked[f >> 1]) return;
    marked[f >> 1] = 1;
    mark_live(nodes[f >> 1].lo);
    mark_live(nodes[f >> 1].hi);
}

static Ref remap_live(Ref f) { return f < 2 ? f : (compact_map[f >> 1] << 1) | (f & 1); }

static void collect_live(Ref *excluded, Ref *losing, Ref *winning) {
    if (!gc_enabled || node_count < node_limit/8) return;
    memset(marked, 0, node_count);
    mark_live(domain); mark_live(*excluded); mark_live(*losing); mark_live(*winning);
    for (unsigned v = 0; v < variables; v++) mark_live(xvars[v]);
    for (unsigned k = 0; k < certificate_count; k++) mark_live(certificate_sets[k]);
#ifdef STRIP3_PARTIAL_KERNEL_SYMBOLIC
    for (unsigned k = 0; k < certificate_count; k++) mark_live(partial_winning[k]);
#endif
    size_t next = 1;
    /* Children are older than parents; in-place compaction is safe. */
    for (size_t id = 1; id < node_count; id++) if (marked[id]) {
        struct node node = nodes[id];
        assert(node.lo < 2 || (node.lo >> 1) < id);
        assert(node.hi < 2 || (node.hi >> 1) < id);
        compact_map[id] = (uint32_t)next;
        nodes[next++] = (struct node){remap_live(node.lo), remap_live(node.hi), node.var};
    }
    domain = remap_live(domain); *excluded = remap_live(*excluded);
    *losing = remap_live(*losing); *winning = remap_live(*winning);
    for (unsigned v = 0; v < variables; v++) xvars[v] = remap_live(xvars[v]);
    for (unsigned k = 0; k < certificate_count; k++) certificate_sets[k] = remap_live(certificate_sets[k]);
#ifdef STRIP3_PARTIAL_KERNEL_SYMBOLIC
    for (unsigned k = 0; k < certificate_count; k++) partial_winning[k] = remap_live(partial_winning[k]);
#endif
    node_count = next;
    memset(unique, 0, unique_cap*sizeof *unique);
    memset(cache, 0, cache_cap*sizeof *cache);
    memset(restrict_generation, 0, node_limit*sizeof *restrict_generation);
    generation = 0;
    for (uint32_t id = 1; id < node_count; id++) {
        struct node node = nodes[id];
        size_t index = mix(((uint64_t)node.lo << 32) ^ node.hi ^
                           ((uint64_t)node.var*0x9e3779b97f4a7c15ULL)) & (unique_cap-1);
        while (unique[index]) index = (index+1) & (unique_cap-1);
        unique[index] = id;
    }
    collections++;
}

/* Image of every independent history under S -> V\N[S].  Histories use
 * vertically independent column masks.  Availability in a column is known
 * after its right history neighbour is selected.  This is an exact regular
 * image, not the larger family of arbitrary induced subgraphs. */
static Ref history_memo[101][8][8];
static unsigned char history_seen[101][8][8];
static unsigned strip_width;

static Ref history_image(unsigned column, unsigned left, unsigned here) {
    if (history_seen[column][left][here]) return history_memo[column][left][here];
    history_seen[column][left][here] = 1;
    static const unsigned choices[] = {0, 1, 2, 4, 5};
    Ref result = 0;
    unsigned blocked = here | ((here << 1) & 7) | (here >> 1) | left;
    for (unsigned i = 0; i < 5; i++) {
        unsigned right = choices[i];
        if ((column+1 == strip_width && right) || (right & here)) continue;
        unsigned available = 7 & ~(blocked | right);
        Ref term = column+1 == strip_width ? 1 : history_image(column+1, here, right);
        /* Descending variable order prepends one output column to its suffix. */
        for (int r = 2; r >= 0; r--) {
            unsigned var = 3*column+(unsigned)r;
            term = (available >> r) & 1 ? make_node(var, 0, term) : make_node(var, term, 0);
        }
        result = bdd_or(result, term);
    }
    history_memo[column][left][here] = result;
    return result;
}

static Ref reachable_domain(unsigned n) {
    static const unsigned choices[] = {0, 1, 2, 4, 5};
    strip_width = n;
    Ref result = 0;
    for (unsigned i = 0; i < 5; i++) result = bdd_or(result, history_image(0, 0, choices[i]));
    return result;
}

/* The checker only needs the final kernels, not any iteration history.
 * This producer completes its lower kernels; the certificate semantics also
 * permit partial lower kernels with locally checked Nim-move witnesses. */
static int write_kernel_certificate(unsigned n, unsigned k, Ref last) {
    if (!certificate_path) return 1;
    if (!compact_map) compact_map = calloc(node_limit, sizeof *compact_map);
    if (!marked) marked = calloc(node_limit, sizeof *marked);
    if (!compact_map || !marked) return 0;
    memset(marked, 0, node_count);
    mark_live(domain); mark_live(last);
    for (unsigned j = 0; j < k; j++) mark_live(certificate_sets[j]);
    size_t next = 1;
    for (size_t id = 1; id < node_count; id++) if (marked[id]) compact_map[id] = (uint32_t)next++;
    FILE *file = fopen(certificate_path, "w");
    if (!file) return 0;
    fprintf(file, "node-kayles-symbolic-kernel-v1\n%u %u %zu\n%u\n", n, k, next, remap_live(domain));
    for (unsigned j = 0; j <= k; j++) fprintf(file, "%u%c", remap_live(j == k ? last : certificate_sets[j]), j == k ? '\n' : ' ');
    for (size_t id = 1; id < node_count; id++) if (marked[id]) {
        struct node node = nodes[id];
        fprintf(file, "%u %u %u\n", node.var, remap_live(node.lo), remap_live(node.hi));
    }
    int okay = !ferror(file);
    if (fclose(file)) okay = 0;
    fprintf(stderr, "kernel_certificate nodes=%zu roots=%u\n", next-1, k+1);
    return okay;
}
#endif

int main(int argc, char **argv) {
    int domain_only = 0;
#ifdef STRIP3_REACHABLE_SYMBOLIC
    if (argc == 5 && !strcmp(argv[1], "board") && !strcmp(argv[3], "--certificate")) {
        certificate_path = argv[4]; argc = 3;
    }
    domain_only = argc == 3 && !strcmp(argv[1], "domain");
#endif
    if (argc != 3 || (!domain_only && strcmp(argv[1], "board") && strcmp(argv[1], "table"))) return 2;
    char *end;
    long n = strtol(argv[2], &end, 10);
    int table = !strcmp(argv[1], "table");
    if (!*argv[2] || *end || n < 1 || n > 100 || ((table || domain_only) && n > 6)) return 2;
    variables = (unsigned)(3*n);
    const char *env = getenv("STRIP3_SYMBOLIC_BITS");
    long bits = env ? strtol(env, &end, 10) : 20;
    if (bits < 8 || bits > 23 || (env && (!*env || *end))) return 2;
    node_limit = (size_t)1 << bits; unique_cap = node_limit*2; cache_cap = node_limit;
    nodes = calloc(node_limit, sizeof *nodes);
    unique = calloc(unique_cap, sizeof *unique);
    cache = calloc(cache_cap, sizeof *cache);
    restrict_generation = calloc(node_limit, sizeof *restrict_generation);
    restrict_cache = calloc(node_limit, sizeof *restrict_cache);
    if (!nodes || !unique || !cache || !restrict_generation || !restrict_cache) return 3;
    int *values = NULL;
#ifdef STRIP3_REACHABLE_SYMBOLIC
    const char *gc_env = getenv("STRIP3_SYMBOLIC_GC");
    if (gc_env && strcmp(gc_env, "0") && strcmp(gc_env, "1")) return 2;
    gc_enabled = !gc_env || !strcmp(gc_env, "1");
    if (gc_enabled) {
        compact_map = calloc(node_limit, sizeof *compact_map);
        marked = calloc(node_limit, sizeof *marked);
        if (!compact_map || !marked) return 3;
    }
#endif
    for (unsigned v = 0; v < variables; v++) xvars[v] = make_node(v, 0, 1);
#ifdef STRIP3_REACHABLE_SYMBOLIC
    domain = reachable_domain((unsigned)n);
    fprintf(stderr, "reachable_domain nodes=%zu apply_calls=%llu\n", node_count,
            (unsigned long long)apply_calls);
    if (domain_only) {
        for (uint64_t m = 0; m < ((uint64_t)1 << variables); m++) printf("%d\n", evaluate(domain, m, 0));
        goto done;
    }
#endif
    size_t masks = table ? (size_t)1 << variables : 0, assigned = 0, wanted = 0;
    if (table) {
        values = malloc(masks*sizeof *values);
        if (!values) return 3;
        for (size_t m = 0; m < masks; m++) {
            values[m] = -1;
            wanted += evaluate(domain, m, 0);
        }
    }
    Ref excluded = 0;
    int g;
#ifdef STRIP3_PARTIAL_KERNEL_SYMBOLIC
    const char *initial_env = getenv("STRIP3_PARTIAL_INITIAL");
    long initial = initial_env ? strtol(initial_env, &end, 10) : 7;
    if (initial < 0 || initial > 300 || (initial_env && (!*initial_env || *end))) return 2;
    unsigned ceiling = (unsigned)initial < variables ? (unsigned)initial : variables;
    certificate_count = ceiling+1;
    for (unsigned round = 1; ; round++) {
        Ref dummy_a = 0, dummy_b = 0, dummy_c = 0;
        collect_live(&dummy_a, &dummy_b, &dummy_c);
        int changed = 0;
        Ref prefix = 0;
        for (unsigned k = 0; k <= ceiling; k++) {
            Ref winning = bdd_or(prefix, predecessor(certificate_sets[k]));
            changed |= winning != partial_winning[k];
            partial_winning[k] = winning;
            prefix = bdd_or(prefix, certificate_sets[k]);
        }
        Ref smaller_false = domain;
        for (unsigned k = 0; k <= ceiling; k++) {
            Ref losing = bdd_and(smaller_false, predecessor(partial_winning[k]^1)^1);
            changed |= losing != certificate_sets[k];
            certificate_sets[k] = losing;
            smaller_false = bdd_and(smaller_false, partial_winning[k]);
            if (!table && evaluate(losing, 0, 1)) {
                fprintf(stderr, "partial round=%u ceiling=%u k=%u nodes=%zu apply_calls=%llu\n",
                        round, ceiling, k, node_count, (unsigned long long)apply_calls);
                if (!write_kernel_certificate((unsigned)n, k, losing)) {
                    fprintf(stderr, "symbolic: certificate write failed; incomplete\n"); return 3;
                }
                printf("3\t%ld\t%u\n", n, k);
                goto done;
            }
        }
        fprintf(stderr, "partial round=%u ceiling=%u nodes=%zu apply_calls=%llu\n",
                round, ceiling, node_count, (unsigned long long)apply_calls);
        if (table) {
            for (size_t m = 0; m < masks; m++) if (values[m] < 0 && evaluate(domain, m, 0)) {
                for (unsigned k = 0; k <= ceiling; k++) if (evaluate(certificate_sets[k], m, 0)) {
                    values[m] = (int)k; assigned++; break;
                }
            }
            if (assigned == wanted) {
                for (size_t m = 0; m < masks; m++) printf("%d\n", values[m]);
                goto done;
            }
        }
        /* No assumption that SG<=the initial ceiling: expand on a completed
         * refutation of the current board range, or on a saturated table. */
        int expand = table ? !changed : evaluate(smaller_false, 0, 1);
        if (expand && ceiling < variables) {
            ceiling++;
            certificate_count = ceiling+1;
        } else if (!changed) {
            fprintf(stderr, "partial: invalid saturated classification; incomplete\n"); return 1;
        }
    }
#endif
    for (g = 0; g <= (int)variables; g++) {
        Ref losing = 0, winning = excluded;
        unsigned steps = 0;
        do {
#ifdef STRIP3_REACHABLE_SYMBOLIC
            collect_live(&excluded, &losing, &winning);
#endif
            Ref unexcluded = domain == 1 ? excluded^1 : bdd_and(domain, excluded^1);
            losing = bdd_and(unexcluded, predecessor(winning^1)^1);
            winning = bdd_or(excluded, predecessor(losing));
            assert(!bdd_and(losing, winning));
            steps++;
            /* Every graph move reduces order; at most |V|+1 layers. */
            if (steps > variables+1) { fprintf(stderr, "symbolic: invalid fixed point\n"); return 1; }
            if (!table && evaluate(losing, 0, 1)) {
                fprintf(stderr, "symbolic k=%d steps=%u nodes=%zu apply_calls=%llu\n", g, steps, node_count,
                        (unsigned long long)apply_calls);
#ifdef STRIP3_REACHABLE_SYMBOLIC
                if (!write_kernel_certificate((unsigned)n, (unsigned)g, losing)) {
                    fprintf(stderr, "symbolic: certificate write failed; incomplete\n"); return 3;
                }
#endif
                printf("3\t%ld\t%d\n", n, g);
                goto done;
            }
        } while (winning != (domain == 1 ? losing^1 : bdd_and(domain, losing^1)));
        fprintf(stderr, "symbolic k=%d steps=%u nodes=%zu apply_calls=%llu\n", g, steps, node_count,
                (unsigned long long)apply_calls);
        if (table) for (size_t m = 0; m < masks; m++) if (evaluate(losing, m, 0)) {
            assert(values[m] == -1);
            values[m] = g; assigned++;
        }
        excluded = bdd_or(excluded, losing);
#ifdef STRIP3_REACHABLE_SYMBOLIC
        if (certificate_path) certificate_sets[certificate_count++] = losing;
#endif
        if (table && assigned == wanted) break;
    }
    if (table && assigned == wanted) for (size_t m = 0; m < masks; m++) printf("%d\n", values[m]);
    else { fprintf(stderr, "symbolic: incomplete SG classification\n"); return 1; }
done:
#ifdef STRIP3_REACHABLE_SYMBOLIC
    fprintf(stderr, "reachable_stats created_nodes=%llu peak_nodes=%zu retained_nodes=%zu gc=%llu apply_calls=%llu\n",
            (unsigned long long)created_nodes, peak_nodes, node_count, (unsigned long long)collections,
            (unsigned long long)apply_calls);
    free(compact_map); free(marked);
#endif
    free(values); free(nodes); free(unique); free(cache); free(restrict_generation); free(restrict_cache);
    return 0;
}
