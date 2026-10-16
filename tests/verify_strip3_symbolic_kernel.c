/* Standalone symbolic SG-kernel checker: no solver source or solver calls.
 *
 * Uses ordinary (uncomplemented) BDDs. Reconstructs the reachable domain by
 * existentially quantifying an independent-set/availability Boolean relation,
 * rather than importing the producer's column-history transducer.
 * Input: plaintext node-kayles-symbolic-kernel-v1 certificate on stdin.
 */
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint32_t B;
struct vertex { unsigned variable; B low, high, negative; };
struct memo { B a, b, result; unsigned operation; };
static struct vertex *graph;
static struct memo *memo;
static uint32_t *unique, *quantified, *restricted, *restriction_stamp;
static size_t used = 2, capacity, unique_size, memo_size;
static unsigned stamp, width, cells, levels;
static int history_mode;
static unsigned char removed[600];
static B literals[600], reachable;
static uint64_t applications, predecessor_calls;

static void require(int condition, const char *message) {
    if (!condition) { fprintf(stderr, "kernel check failed: %s\n", message); exit(1); }
}

static uint64_t scramble(uint64_t a) {
    a ^= a >> 33; a *= UINT64_C(0xff51afd7ed558ccd);
    a ^= a >> 33; a *= UINT64_C(0xc4ceb9fe1a85ec53);
    return a ^ (a >> 33);
}

static unsigned variable(B a) { return a < 2 ? levels : graph[a].variable; }

static B node(unsigned v, B low, B high) {
    if (low == high) return low;
    require(v < variable(low) && v < variable(high), "unordered internal BDD");
    size_t slot = scramble(((uint64_t)low << 32) ^ high ^ ((uint64_t)v << 17)) & (unique_size-1);
    while (unique[slot]) {
        B a = unique[slot];
        if (graph[a].variable == v && graph[a].low == low && graph[a].high == high) return a;
        slot = (slot+1) & (unique_size-1);
    }
    require(used < capacity, "checker node budget exhausted");
    B a = (B)used++;
    graph[a] = (struct vertex){v, low, high, 0};
    unique[slot] = a;
    return a;
}

static B negate(B a) {
    if (a < 2) return a ^ 1;
    if (graph[a].negative) return graph[a].negative;
    struct vertex copy = graph[a];
    B low = negate(copy.low), high = negate(copy.high);
    B answer = node(copy.variable, low, high);
    graph[a].negative = answer; graph[answer].negative = a;
    return answer;
}

/* operation 1 is intersection, 2 is union. */
static B apply(unsigned operation, B a, B b) {
    applications++;
    if (a == b) return a;
    if (operation == 1) {
        if (!a || !b) return 0;
        if (a == 1) return b;
        if (b == 1) return a;
    } else {
        if (a == 1 || b == 1) return 1;
        if (!a) return b;
        if (!b) return a;
    }
    if (a > b) { B t = a; a = b; b = t; }
    size_t slot = scramble(((uint64_t)a << 32) ^ b ^ ((uint64_t)operation << 23)) & (memo_size-1);
    struct memo entry = memo[slot];
    if (entry.operation == operation && entry.a == a && entry.b == b) return entry.result;
    unsigned va = variable(a), vb = variable(b), v = va < vb ? va : vb;
    B al = a, ah = a, bl = b, bh = b;
    if (va == v) { al = graph[a].low; ah = graph[a].high; }
    if (vb == v) { bl = graph[b].low; bh = graph[b].high; }
    B low = apply(operation, al, bl), high = apply(operation, ah, bh);
    B answer = node(v, low, high);
    memo[slot] = (struct memo){a, b, answer, operation};
    return answer;
}

static B intersect(B a, B b) { return apply(1, a, b); }
static B unite(B a, B b) { return apply(2, a, b); }

static B quantify_history(B a) {
    if (a < 2) return a;
    if (quantified[a]) return quantified[a]-1;
    struct vertex copy = graph[a];
    B low = quantify_history(copy.low), high = quantify_history(copy.high);
    B answer = copy.variable & 1 ? unite(low, high) : node(copy.variable, low, high);
    quantified[a] = answer+1;
    return answer;
}

static int distance(unsigned u, unsigned v) {
    return abs((int)(u/3)-(int)(v/3)) + abs((int)(u%3)-(int)(v%3));
}

static B independent_reachable_domain(void) {
    B relation = 1;
    for (unsigned u = 0; u < cells; u++) for (unsigned v = u+1; v < cells; v++) {
        if (distance(u, v) == 1)
            relation = intersect(relation, unite(negate(literals[2*u+1]), negate(literals[2*v+1])));
    }
    for (unsigned v = 0; v < cells; v++) {
        B blocked = 0;
        for (unsigned u = 0; u < cells; u++) if (distance(u, v) <= 1)
            blocked = unite(blocked, literals[2*u+1]);
        B available = literals[2*v];
        B equation = unite(intersect(available, negate(blocked)), intersect(negate(available), blocked));
        relation = intersect(relation, equation);
    }
    return quantify_history(relation);
}

static B independent_history_domain(void) {
    B answer = 1;
    for (unsigned u = 0; u < cells; u++) for (unsigned v = u+1; v < cells; v++)
        if (distance(u, v) == 1)
            answer = intersect(answer, unite(negate(literals[2*u]), negate(literals[2*v])));
    return answer;
}

static B delete_neighbourhood(B a) {
    if (a < 2) return a;
    if (restriction_stamp[a] == stamp) return restricted[a];
    struct vertex copy = graph[a];
    B answer;
    if (removed[copy.variable])
        answer = delete_neighbourhood(removed[copy.variable] == 2 ? copy.high : copy.low);
    else {
        B low = delete_neighbourhood(copy.low), high = delete_neighbourhood(copy.high);
        answer = node(copy.variable, low, high);
    }
    restriction_stamp[a] = stamp; restricted[a] = answer;
    return answer;
}

static B predecessor(B target) {
    predecessor_calls++;
    B answer = 0;
    for (unsigned v = 0; v < cells; v++) {
        memset(removed, 0, sizeof removed);
        B legal = literals[2*v];
        if (history_mode) {
            legal = 1;
            for (unsigned u = 0; u < cells; u++) if (distance(u, v) <= 1)
                legal = intersect(legal, negate(literals[2*u]));
            removed[2*v] = 2;
        } else {
            for (unsigned u = 0; u < cells; u++) if (distance(u, v) <= 1) removed[2*u] = 1;
        }
        stamp++;
        require(stamp != 0, "restriction stamp exhausted");
        B after = delete_neighbourhood(target);
        B term = intersect(intersect(reachable, legal), after);
        answer = unite(answer, term);
    }
    return answer;
}

static uint32_t integer(void) {
    char token[64], *end;
    require(scanf("%63s", token) == 1, "truncated certificate");
    require(token[0] >= '0' && token[0] <= '9', "invalid unsigned integer");
    errno = 0;
    unsigned long long value = strtoull(token, &end, 10);
    require(!errno && !*end && value <= UINT32_MAX, "integer out of range");
    return (uint32_t)value;
}

static B from_reference(B ref, const B *translated) {
    B answer = translated[ref >> 1];
    return ref & 1 ? negate(answer) : answer;
}

int main(void) {
    char header[80];
    require(fgets(header, sizeof header, stdin) != NULL, "missing format");
    history_mode = !strcmp(header, "node-kayles-history-kernel-v1\n");
    require(history_mode || !strcmp(header, "node-kayles-symbolic-kernel-v1\n"), "wrong format");
    width = integer(); unsigned target = integer(), count = integer();
    require(width >= 1 && width <= 100 && target <= 3*width, "invalid board or nimber");
    cells = 3*width; levels = 2*cells;
    const char *setting = getenv("STRIP3_KERNEL_BITS");
    char *end;
    long bits = setting ? strtol(setting, &end, 10) : 20;
    require(bits >= 8 && bits <= 22 && (!setting || (*setting && !*end)), "invalid checker budget");
    capacity = (size_t)1 << bits; unique_size = capacity*2; memo_size = capacity;
    require(count >= 1 && count <= capacity, "invalid source node count");
    graph = calloc(capacity, sizeof *graph); unique = calloc(unique_size, sizeof *unique);
    memo = calloc(memo_size, sizeof *memo); quantified = calloc(capacity, sizeof *quantified);
    restricted = calloc(capacity, sizeof *restricted); restriction_stamp = calloc(capacity, sizeof *restriction_stamp);
    B *translated = calloc(count, sizeof *translated);
    unsigned *source_levels = calloc(count, sizeof *source_levels);
    require(graph && unique && memo && quantified && restricted && restriction_stamp && translated && source_levels, "allocation failed");
    B domain_ref = integer(), roots[301];
    require((domain_ref >> 1) < count, "invalid domain reference");
    for (unsigned k = 0; k <= target; k++) { roots[k] = integer(); require((roots[k] >> 1) < count, "invalid kernel reference"); }
    source_levels[0] = cells;
    for (unsigned id = 1; id < count; id++) {
        unsigned v = integer(); B low = integer(), high = integer();
        require(v < cells && (low >> 1) < id && (high >> 1) < id && low != high && !(low & 1), "invalid source BDD node");
        require(source_levels[low >> 1] > v && source_levels[high >> 1] > v, "unordered source BDD");
        source_levels[id] = v;
        B a = from_reference(low, translated), b = from_reference(high, translated);
        translated[id] = node(2*v, a, b);
    }
    require(scanf(" %1s", header) == EOF, "trailing certificate data");
    B claimed_domain = from_reference(domain_ref, translated), kernels[301];
    for (unsigned k = 0; k <= target; k++) kernels[k] = from_reference(roots[k], translated);
    for (unsigned v = 0; v < levels; v++) literals[v] = node(v, 0, 1);
    reachable = history_mode ? independent_history_domain() : independent_reachable_domain();
    require(reachable == claimed_domain, "domain differs from its independent reconstruction");
    B excluded = 0, witnesses[301];
    unsigned complete_lower_kernels = 0;
    for (unsigned k = 0; k <= target; k++) {
        B kernel = kernels[k];
        require(!intersect(kernel, negate(reachable)), "kernel outside the reachable domain");
        require(!intersect(kernel, excluded), "kernel overlaps a lower nimber");
        B can_reach = predecessor(kernel), winning = unite(excluded, can_reach);
        witnesses[k] = winning;
        require(!intersect(kernel, winning), "kernel has a move to itself or a lower classification");
        /* A move in the auxiliary Nim heap must also have a winning witness.
         * Lower kernels may be partial: no global classification is assumed. */
        for (unsigned j = 0; j < k; j++)
            require(!intersect(kernel, negate(witnesses[j])), "Nim move lacks a winning witness");
        B unknown = intersect(reachable, negate(winning));
        B escape = unknown == kernel ? can_reach : predecessor(unknown);
        require(!intersect(kernel, escape), "kernel has a successor outside its winning witnesses");
        if (k < target && unknown == kernel) complete_lower_kernels++;
        excluded = unite(excluded, kernel);
    }
    B board = kernels[target];
    while (board >= 2) board = history_mode ? graph[board].low : graph[board].high;
    require(board == 1, "full board is absent from the last kernel");
    printf("{\"n\":%u,\"SG\":%u,\"certificate_nodes\":%u,\"checker_nodes\":%zu,\"complete_lower_kernels\":%u,\"predecessor_calls\":%llu,\"apply_calls\":%llu}\n",
           width, target, count-1, used-2, complete_lower_kernels, (unsigned long long)predecessor_calls, (unsigned long long)applications);
    free(graph); free(unique); free(memo); free(quantified); free(restricted); free(restriction_stamp);
    free(translated); free(source_levels);
    return 0;
}
