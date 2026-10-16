/* OFFLINE proof diagnostic; full SG oracle cost is part of construction.
 * Minimise an unfolded equality proof TREE in a stated calculus, or use the
 * first oracle-confirmed winning branch. Merge identical (component,target)
 * nodes afterwards. This does NOT minimise the shared DAG and is NOT an
 * algorithm for unknown large strips. Exported proofs need no SG oracle.
 * Rules: path/cycle/safe-rotation axioms; every graph and smaller-Nim option
 * for equality; one winning graph/Nim option for inequality. Disjoint sums
 * cancel equal components and certify every nonlargest component exactly.
 */
#ifndef STRIP3_PROOF_EXTERNAL_BASE
#define main oracle_cli_main
#include "strip3.c"
#undef main
#endif

enum { R_BOUND=1, R_PATH, R_CYCLE, R_ROTATION, R_ALL, R_MOVE, R_NIM };
struct proof_entry {
    Mask key, cost;
    uint32_t id;
    short move;
    unsigned char k, rule, truth;
};
struct option { Mask mask; uint64_t signature; short move; unsigned char g; };
static struct proof_entry *proofs;
static size_t proof_cap, proof_count, selected_count, selected_cap;
static size_t *selected;
static uint64_t considered_moves, sum_evaluations, emitted_sum_rules, emitted_edges;
static int first_policy;
static int (*proof_value_hook)(Mask);
static int (*proof_truth_hook)(Mask, int);
static const Mask COST_INF = (Mask)-1;

static const char *policy_name(void) {
    return proof_value_hook ? (first_policy ? "available_first" : "available_tree_min")
                            : (first_policy ? "first" : "tree_min");
}

static int proof_value(Mask m) { return proof_value_hook ? proof_value_hook(m) : component_value(m); }

static int proof_truth(Mask m, int k) {
    return proof_truth_hook ? proof_truth_hook(m, k) : proof_value(m) == k;
}

static Mask add_cost(Mask a, Mask b) { return COST_INF-a < b ? COST_INF : a+b; }

static void decimal(FILE *out, Mask x) {
    char text[40]; int n = 0;
    do { text[n++] = (char)('0'+x%10); x /= 10; } while (x);
    while (n) fputc(text[--n], out);
}

static void hexmask(FILE *out, Mask m) {
    if (m >> 64) fprintf(out, "%llx%016llx", (unsigned long long)(m >> 64), (unsigned long long)m);
    else fprintf(out, "%llx", (unsigned long long)m);
}

static size_t proof_slot(Mask m, int k) {
    size_t i = (hash(m)^((uint64_t)k*0x9e3779b97f4a7c15ULL)) & (proof_cap-1);
    while (proofs[i].key && (proofs[i].key != m || proofs[i].k != k)) i = (i+1)&(proof_cap-1);
    if (!proofs[i].key) {
        if (proof_count >= proof_cap*7/10) {
            fprintf(stderr, "proof-cost: proof budget reached; incomplete\n"); exit(3);
        }
        proofs[i].key = m; proofs[i].k = (unsigned char)k; proof_count++;
    }
    return i;
}

static int parts_of(Mask m, Mask *parts) {
    int count = 0;
    while (m) {
        Mask part = (Mask)1 << lowbit(m), frontier = part;
        while (frontier) { frontier = neighbours(frontier)&m&~part; part |= frontier; }
        m &= ~part;
        Mask key = canonical(part); int found = -1;
        for (int i = 0; i < count; i++) if (parts[i] == key) { found = i; break; }
        if (found >= 0) parts[found] = parts[--count];
        else parts[count++] = key;
    }
    /* Largest first; ties by increasing mask. This fixed decomposition rule
     * is part of the proof calculus whose minimum TREE is being measured. */
    for (int i = 1; i < count; i++) {
        Mask x = parts[i]; int j = i;
        while (j && (order(parts[j-1]) < order(x) ||
               (order(parts[j-1]) == order(x) && parts[j-1] > x))) {
            parts[j] = parts[j-1]; j--;
        }
        parts[j] = x;
    }
    return count;
}

static uint64_t tuple_hash(Mask *parts, int count) {
    uint64_t h = (uint64_t)count;
    for (int i = 0; i < count; i++) h = (h^hash(parts[i]))*0xbf58476d1ce4e5b9ULL;
    return h;
}

static int same_sum(Mask m, Mask *parts, int count) {
    Mask other[MAXCELLS]; int nc = parts_of(m, other);
    return nc == count && !memcmp(parts, other, count*sizeof *parts);
}

static int option_list(Mask m, struct option *options) {
    int count = 0, n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    for (int pass = 0; pass < 2; pass++) for (int step = 0; step < n; step++) {
        int c = n/2+(step%2 ? -(step+1)/2 : step/2);
        for (int r = 0; r < h; r++) {
            if ((r == 1) != (pass == 0)) continue;
            int v = 3*c+r; Mask bit = (Mask)1 << v;
            if (!(m & bit)) continue;
            Mask child = m & ~(bit | neighbours(bit)), parts[MAXCELLS];
            int np = parts_of(child, parts), duplicate = 0;
            uint64_t sig = tuple_hash(parts, np);
            for (int i = 0; i < count; i++) if (options[i].signature == sig && same_sum(options[i].mask, parts, np)) {
                duplicate = 1; break;
            }
            if (!duplicate) {
                int g = 0;
                for (int i = 0; i < np; i++) {
                    int part_g = proof_value(parts[i]);
                    if (part_g < 0) { g = -1; break; }
                    g ^= part_g;
                }
                options[count++] = (struct option){child, sig, (short)v, (unsigned char)g};
                considered_moves++;
            }
        }
    }
    return count;
}

static int axiom(Mask m) {
    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask both = reverse_bits(m, 3*n) >> (3-h);
    if (m == both) {
        int safe = 1;
        for (Mask x = m; x && safe; x &= x-1) {
            int v = lowbit(x), partner = 3*(n-1-v/3)+h-1-v%3;
            if (partner == v || (neighbours((Mask)1 << v) & ((Mask)1 << partner))) safe = 0;
        }
        if (safe) return R_ROTATION;
    }
    int maxdegree = 0, endpoints = 0;
    for (Mask x = m; x; x &= x-1) {
        int degree = order(neighbours((Mask)1 << lowbit(x)) & m);
        if (degree > maxdegree) maxdegree = degree;
        if (degree < 2) endpoints++;
    }
    return maxdegree <= 2 ? (endpoints || order(m) == 1 ? R_PATH : R_CYCLE) : 0;
}

static Mask tree_cost(Mask m, int k);

/* Full-oracle experiments fix pivot 0. With partial search knowledge there
 * can be one component of unknown exact value; choose it as pivot and retain
 * exact proofs for every other component. No missing value is computed here. */
static int sum_plan(Mask m, int k, Mask *parts, int *values, int *pivot, int *target) {
    int count = parts_of(m, parts), unknown = 0, offset = 0;
    *pivot = 0;
    for (int i = 0; i < count; i++) {
        values[i] = proof_value(parts[i]);
        if (values[i] < 0) { unknown++; *pivot = i; }
    }
    if (unknown > 1) return -1;
    for (int i = 0; i < count; i++) if (i != *pivot) offset ^= values[i];
    *target = k^offset;
    if (count && proof_truth(parts[*pivot], *target) < 0) return -1;
    return count;
}

static Mask sum_cost(Mask m, int k) {
    sum_evaluations++;
    Mask parts[MAXCELLS]; int values[MAXCELLS], pivot, target;
    int count = sum_plan(m, k, parts, values, &pivot, &target);
    if (count < 0) return COST_INF;
    Mask cost = 1;
    for (int i = 0; i < count; i++) cost = add_cost(cost, tree_cost(parts[i], i == pivot ? target : values[i]));
    return cost;
}

static Mask tree_cost(Mask m, int k) {
    size_t slot = proof_slot(m, k);
    struct proof_entry *e = &proofs[slot];
    if (e->cost) return e->cost;
    int g = proof_value(m), nodes = order(m), truth = proof_truth(m, k);
    if (truth < 0) { e->cost = COST_INF; return COST_INF; }
    e->truth = (unsigned char)truth;
    int primitive = k > nodes ? R_BOUND : axiom(m);
    if (primitive) { e->rule = (unsigned char)primitive; e->cost = 1; return 1; }
    struct option options[MAXCELLS]; int count = option_list(m, options);
    Mask cost = COST_INF;
    if (truth) {
        cost = 1; e->rule = R_ALL;
        for (int i = 0; i < count; i++) cost = add_cost(cost, sum_cost(options[i].mask, k));
        for (int j = 0; j < k; j++) cost = add_cost(cost, tree_cost(m, j));
    } else {
        for (int i = 0; i < count; i++) if (options[i].g == k) {
            Mask candidate = add_cost(1, sum_cost(options[i].mask, k));
            if (candidate < cost) { cost = candidate; e->rule = R_MOVE; e->move = options[i].move; }
            if (first_policy) break;
        }
        if (g >= 0 && g < k && (!first_policy || cost == COST_INF)) {
            Mask candidate = add_cost(1, tree_cost(m, g));
            if (candidate < cost) { cost = candidate; e->rule = R_NIM; e->move = (short)g; }
        }
    }
    e->cost = cost;
    if (!e->rule && !proof_value_hook) { fprintf(stderr, "proof-cost: no finite proof; incomplete\n"); exit(3); }
    return cost;
}

static uint32_t select_query(Mask m, int k);

static void select_sum(Mask m, int k) {
    Mask parts[MAXCELLS]; int values[MAXCELLS], pivot, target;
    int count = sum_plan(m, k, parts, values, &pivot, &target);
    if (count < 0) { fprintf(stderr, "proof-cost: unavailable sum proof; incomplete\n"); exit(3); }
    for (int i = 0; i < count; i++) select_query(parts[i], i == pivot ? target : values[i]);
}

static uint32_t select_query(Mask m, int k) {
    size_t slot = proof_slot(m, k); struct proof_entry *e = &proofs[slot];
    if (e->id) return e->id;
    if (!e->cost) tree_cost(m, k);
    if (selected_count == selected_cap) {
        size_t next = selected_cap ? selected_cap*2 : 1024;
        size_t *p = realloc(selected, next*sizeof *p);
        if (!p) { fprintf(stderr, "proof-cost: certificate allocation failed; incomplete\n"); exit(3); }
        selected = p; selected_cap = next;
    }
    e->id = (uint32_t)(selected_count+1); selected[selected_count++] = slot;
    if (e->rule == R_ALL) {
        struct option options[MAXCELLS]; int count = option_list(m, options);
        for (int i = 0; i < count; i++) select_sum(options[i].mask, k);
        for (int j = 0; j < k; j++) select_query(m, j);
    } else if (e->rule == R_MOVE) {
        Mask bit = (Mask)1 << e->move; select_sum(m & ~(bit | neighbours(bit)), k);
    } else if (e->rule == R_NIM) select_query(m, e->move);
    return e->id;
}

static void emit_sum(FILE *out, Mask m, int k) {
    Mask parts[MAXCELLS]; int values[MAXCELLS], pivot, target;
    int count = sum_plan(m, k, parts, values, &pivot, &target);
    if (count < 0) { fprintf(stderr, "proof-cost: unavailable emitted sum; incomplete\n"); exit(3); }
    int truth = count ? proof_truth(parts[pivot], target) : k == 0;
    fprintf(out, "{\"k\":%d,\"equal\":%s,\"pivot\":%d,\"parts\":[", k, truth ? "true" : "false", pivot);
    for (int i = 0; i < count; i++) { if (i) fputc(',', out); fputc('"', out); hexmask(out, parts[i]); fputc('"', out); }
    fprintf(out, "],\"children\":[");
    for (int i = 0; i < count; i++) {
        int child_target = i == pivot ? target : values[i];
        if (i) fputc(',', out);
        fprintf(out, "%u", proofs[proof_slot(parts[i], child_target)].id);
        emitted_edges++;
    }
    fprintf(out, "]}"); emitted_sum_rules++;
}

static void emit_certificate(const char *path, int n, Mask root, int g, Mask cost) {
    uint32_t root_id = select_query(root, g);
    FILE *out = fopen(path, "w");
    if (!out) { perror("certificate"); exit(2); }
    fprintf(out, "{\"format\":\"strip3-equality-proof-v1\",\"policy\":\"%s\",\"n\":%d,\"SG\":%d,\"root\":%u,\"tree_cost\":\"",
            policy_name(), n, g, root_id);
    decimal(out, cost); fprintf(out, "\",\"nodes\":[\n");
    for (size_t i = 0; i < selected_count; i++) {
        struct proof_entry *e = &proofs[selected[i]]; Mask m = e->key; int k = e->k;
        if (i) fprintf(out, ",\n");
        fprintf(out, "{\"id\":%u,\"mask\":\"", e->id); hexmask(out, m);
        fprintf(out, "\",\"k\":%d,\"equal\":%s,\"cost\":\"", k, e->truth ? "true" : "false");
        decimal(out, e->cost); fprintf(out, "\",\"rule\":%d", e->rule);
        if (e->rule == R_ALL) {
            struct option options[MAXCELLS]; int count = option_list(m, options);
            fprintf(out, ",\"options\":[");
            for (int j = 0; j < count; j++) { if (j) fputc(',', out); emit_sum(out, options[j].mask, k); }
            fprintf(out, "],\"lower\":[");
            for (int j = 0; j < k; j++) { if (j) fputc(',', out); fprintf(out, "%u", proofs[proof_slot(m, j)].id); emitted_edges++; }
            fprintf(out, "]");
        } else if (e->rule == R_MOVE) {
            Mask bit = (Mask)1 << e->move;
            fprintf(out, ",\"move\":%d,\"witness\":", e->move); emit_sum(out, m & ~(bit | neighbours(bit)), k);
        } else if (e->rule == R_NIM) {
            fprintf(out, ",\"lower_k\":%d,\"child\":%u", e->move, proofs[proof_slot(m, e->move)].id); emitted_edges++;
        }
        fprintf(out, "}");
    }
    fprintf(out, "\n]}\n");
    if (fclose(out)) { perror("certificate close"); exit(2); }
}

static int integer(const char *text, int lo, int hi) {
    char *end; long value = strtol(text, &end, 10);
    return *text && !*end && value >= lo && value <= hi ? (int)value : -1;
}

#ifndef STRIP3_PROOF_EXTERNAL_MAIN
int main(int argc, char **argv) {
    int n = 0, batch = argc == 2 && !strcmp(argv[1], "queries");
    if (!batch && !(argc == 5 && !strcmp(argv[1], "board") && !strcmp(argv[3], "--certificate"))) goto usage;
    if (!batch && (n = integer(argv[2], 1, MAXCOLS)) < 0) goto usage;
    const char *env = getenv("STRIP3_PROOF_FIRST");
    if (env && (first_policy = integer(env, 0, 1)) < 0) goto usage;
    env = getenv("STRIP3_PROOF_BITS"); int bits = env ? integer(env, 8, 24) : 20;
    if (bits < 0) goto usage;
    proof_cap = (size_t)1 << bits; proofs = calloc(proof_cap, sizeof *proofs);
    if (!proofs) { fprintf(stderr, "proof-cost: allocation failed; incomplete\n"); return 3; }
    /* Deliberately avoid finite comb tables in the diagnostic calculus. */
    setenv("STRIP3_THREADS", "1", 1); setenv("STRIP3_REDUCTIONS", "7", 1);
    initialise();
    if (batch) {
        char line[160], hex[40], extra; int k;
        while (fgets(line, sizeof line, stdin)) {
            if (sscanf(line, "%d %39s %d %c", &n, hex, &k, &extra) != 3 || n < 1 || n > MAXCOLS || k < 0 || k > 127 || strlen(hex) > 32) goto usage;
            Mask m = 0;
            for (const char *p = hex; *p; p++) {
                int digit = *p >= '0' && *p <= '9' ? *p-'0' : *p >= 'a' && *p <= 'f' ? *p-'a'+10 : -1;
                if (digit < 0) goto usage;
                m = (m << 4) | digit;
            }
            if (m & ~full_board(n)) goto usage;
            int g = value(m); Mask cost = sum_cost(m, k);
            if (cost == COST_INF) { fprintf(stderr, "proof-cost: tree cost overflow; incomplete\n"); return 3; }
            fprintf(stdout, "{\"equal\":%s,\"g\":%d,\"cost\":\"", g == k ? "true" : "false", g);
            decimal(stdout, cost); fprintf(stdout, "\"}\n");
        }
    } else {
        Mask root = full_board(n); int g = value(root);
        double oracle_time = seconds()-started;
        size_t oracle_states = entries; uint64_t oracle_moves = tst.moves;
        fprintf(stderr, "oracle n=%d SG=%d time=%.6f states=%zu moves=%llu\n", n, g, oracle_time, oracle_states, (unsigned long long)oracle_moves); fflush(stderr);
        Mask cost = tree_cost(root, g);
        if (cost == COST_INF) { fprintf(stderr, "proof-cost: tree cost overflow; incomplete\n"); return 3; }
        emit_certificate(argv[4], n, root, g, cost);
        fprintf(stdout, "{\"n\":%d,\"SG\":%d,\"policy\":\"%s\",\"tree_cost\":\"", n, g, policy_name());
        decimal(stdout, cost);
        fprintf(stdout, "\",\"proof_states\":%zu,\"certificate_queries\":%zu,\"sum_rules\":%llu,\"edges\":%llu,\"considered_moves\":%llu,\"sum_evaluations\":%llu,\"oracle_states\":%zu,\"oracle_moves\":%llu,\"oracle_states_final\":%zu,\"oracle_moves_final\":%llu,\"oracle_calls_final\":%llu,\"oracle_s\":%.6f,\"total_s\":%.6f}\n",
                proof_count, selected_count, (unsigned long long)emitted_sum_rules, (unsigned long long)emitted_edges,
                (unsigned long long)considered_moves, (unsigned long long)sum_evaluations, oracle_states,
                (unsigned long long)oracle_moves, (size_t)entries, (unsigned long long)tst.moves,
                (unsigned long long)tst.calls, oracle_time, seconds()-started);
    }
    return 0;
usage:
    fprintf(stderr, "usage: strip3_proof_cost board N --certificate PATH | queries (N HEX K)\n"); return 2;
}
#endif
