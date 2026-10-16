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
        exit(3);
    }
    uint32_t id = (uint32_t)node_count++;
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
        result = bdd_or(result, bdd_and(xvars[v], restrict_zero(f)));
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

int main(int argc, char **argv) {
    if (argc != 3 || (strcmp(argv[1], "board") && strcmp(argv[1], "table"))) return 2;
    char *end;
    long n = strtol(argv[2], &end, 10);
    int table = !strcmp(argv[1], "table");
    if (!*argv[2] || *end || n < 1 || n > 100 || (table && n > 6)) return 2;
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
    for (unsigned v = 0; v < variables; v++) xvars[v] = make_node(v, 0, 1);
    int *values = NULL;
    size_t masks = table ? (size_t)1 << variables : 0, assigned = 0;
    if (table) {
        values = malloc(masks*sizeof *values);
        if (!values) return 3;
        for (size_t m = 0; m < masks; m++) values[m] = -1;
    }
    Ref excluded = 0;
    int g;
    for (g = 0; g <= (int)variables; g++) {
        Ref losing = 0, winning = excluded;
        unsigned steps = 0;
        do {
            losing = bdd_and(excluded^1, predecessor(winning^1)^1);
            winning = bdd_or(excluded, predecessor(losing));
            assert(!bdd_and(losing, winning));
            steps++;
            /* Every graph move reduces order; at most |V|+1 layers. */
            if (steps > variables+1) { fprintf(stderr, "symbolic: invalid fixed point\n"); return 1; }
            if (!table && evaluate(losing, 0, 1)) {
                fprintf(stderr, "symbolic k=%d steps=%u nodes=%zu apply_calls=%llu\n", g, steps, node_count,
                        (unsigned long long)apply_calls);
                printf("3\t%ld\t%d\n", n, g);
                goto done;
            }
        } while (winning != (losing^1));
        fprintf(stderr, "symbolic k=%d steps=%u nodes=%zu apply_calls=%llu\n", g, steps, node_count,
                (unsigned long long)apply_calls);
        if (table) for (size_t m = 0; m < masks; m++) if (evaluate(losing, m, 0)) {
            assert(values[m] == -1);
            values[m] = g; assigned++;
        }
        excluded = bdd_or(excluded, losing);
        if (table && assigned == masks) break;
    }
    if (table && assigned == masks) for (size_t m = 0; m < masks; m++) printf("%d\n", values[m]);
    else { fprintf(stderr, "symbolic: incomplete SG classification\n"); return 1; }
done:
    free(values); free(nodes); free(unique); free(cache); free(restrict_generation); free(restrict_cache);
    return 0;
}
