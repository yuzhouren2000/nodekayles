/* Depth-first proof-number search for (a disjoint sum)+*k.
 * Node proposition: the position is N. Proof=min(child disproof),
 * disproof=sum(child proof). A solved P singleton proves its exact nimber.
 * Exact canonical tuple keys are compared after hashing. Hard node/key/edge
 * budgets return incomplete rather than forgetting unresolved obligations.
 */
#define TARGET_ENTRY_POINT target_cli_main
#include "strip3_target.c"

enum { INF = 1000000000U };
struct proof_node {
    Mask *parts;
    uint32_t *children;
    uint64_t hash;
    uint32_t proof, disproof;
    unsigned short count, child_count;
    unsigned char k, expanded;
};
static struct proof_node *proof_table;
static Mask *proof_keys;
static uint32_t *proof_edges;
static size_t proof_capacity, proof_entries, key_capacity, key_used, edge_capacity, edge_used;
static uint64_t proof_visits, proof_expansions, proof_generated_moves;
static int small_threshold = 10;

static int proof_single(Mask m, int k);

static void proof_budget(void) {
    fprintf(stderr, "dfpn: budget reached; incomplete; proof_memo=%zu expansions=%llu visits=%llu edges=%zu\n",
            proof_entries, (unsigned long long)proof_expansions, (unsigned long long)proof_visits, edge_used);
    exit(3);
}

static int quick_value(Mask m) {
    struct query_entry *entry = query_slot(m);
    if (entry->exact >= 0) return entry->exact;
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
    int maxdegree = 0, endpoints = 0, nodes = order(m);
    for (Mask x = m; x; x &= x-1) {
        int degree = order(neighbours((Mask)1 << lowbit(x)) & m);
        if (degree > maxdegree) maxdegree = degree;
        if (degree < 2) endpoints++;
        if (maxdegree > 2) break;
    }
    if (maxdegree <= 2) {
        int g = endpoints || nodes == 1 ? path_g[nodes] : (path_g[nodes-3] == 0);
        set_exact(entry, g);
        return g;
    }
    if (nodes <= small_threshold) {
        /* Finite small SG oracle, independently checked by subset mex. */
        component_equality_hook = NULL;
        int (*saved_sum_hook)(Mask *, int, int) = sum_equality_hook;
        sum_equality_hook = NULL;
        int g = exact_value(m);
        sum_equality_hook = saved_sum_hook;
        component_equality_hook = proof_single;
        return g;
    }
    return -1;
}

static void toggle(Mask *parts, int *count, Mask m) {
    for (int i = 0; i < *count; i++) if (parts[i] == m) {
        parts[i] = parts[--*count]; return;
    }
    parts[(*count)++] = m;
}

/* 0 = solved N; 1 = solved P; actual memo slots use indexes >=2. */
static uint32_t proof_key(Mask *input, int count, int k) {
    Mask parts[MAXCELLS]; int kept = 0, total_nodes = 0;
    for (int i = 0; i < count; i++) {
        int g = quick_value(input[i]);
        if (g >= 0) k ^= g;
        else { parts[kept++] = input[i]; total_nodes += order(input[i]); }
    }
    count = kept;
    if (!count) return k == 0 ? 1 : 0;
    if (k > total_nodes) return 0;
    if (count == 1) {
        struct query_entry *entry = query_slot(parts[0]);
        if ((entry->excluded >> k) & 1) return 0;
    }
    for (int i = 1; i < count; i++) {
        Mask x = parts[i]; int j = i;
        while (j && parts[j-1] > x) { parts[j] = parts[j-1]; j--; }
        parts[j] = x;
    }
    uint64_t h = (uint64_t)(count*128+k)*0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < count; i++) h = (h ^ hash(parts[i]))*0xbf58476d1ce4e5b9ULL;
    size_t slot = h & (proof_capacity-1);
    if (slot < 2) slot = 2;
    while (proof_table[slot].count) {
        struct proof_node *e = &proof_table[slot];
        if (e->hash == h && e->count == count && e->k == k && !memcmp(e->parts, parts, count*sizeof *parts))
            return (uint32_t)slot;
        slot = (slot+1) & (proof_capacity-1);
        if (slot < 2) slot = 2;
    }
    if (proof_entries >= proof_capacity*7/10 || key_used+(size_t)count > key_capacity) proof_budget();
    struct proof_node *e = &proof_table[slot];
    e->parts = proof_keys+key_used;
    memcpy(e->parts, parts, count*sizeof *parts);
    key_used += count;
    e->hash = h; e->count = (unsigned short)count; e->k = (unsigned char)k;
    e->proof = e->disproof = 1;
    proof_entries++;
    return (uint32_t)slot;
}

static void append_child(uint32_t *ids, int *count, uint32_t id) {
    for (int i = 0; i < *count; i++) if (ids[i] == id) return;
    ids[(*count)++] = id;
}

static void expand(struct proof_node *e) {
    uint32_t children[MAXCELLS*2+128]; int count = 0;
    for (int a = 0; a < e->count; a++) {
        Mask m = e->parts[a];
        int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
        Mask vertical = flip_rows(m) >> (3-h), both = reverse_bits(m, 3*n) >> (3-h);
        Mask horizontal = flip_rows(both) >> (3-h);
        for (int pass = 0; pass < 2; pass++) for (int step = 0; step < n; step++) {
            int c = n/2 + (step%2 ? -(step+1)/2 : step/2);
            for (int r = 0; r < h; r++) {
                int v = 3*c+r;
                if (!(m >> v & 1) || ((r == 1) != (pass == 0))) continue;
                if (m == vertical && 3*c+h-1-r < v) continue;
                if (m == horizontal && 3*(n-1-c)+r < v) continue;
                if (m == both && 3*(n-1-c)+h-1-r < v) continue;
                Mask next[MAXCELLS], parts[MAXCELLS]; int nc = 0;
                for (int i = 0; i < e->count; i++) if (i != a) toggle(next, &nc, e->parts[i]);
                Mask bit = (Mask)1 << v;
                int cn = split_position(m & ~(bit | neighbours(bit)), parts);
                for (int i = 0; i < cn; i++) toggle(next, &nc, parts[i]);
                uint32_t id = proof_key(next, nc, e->k);
                proof_generated_moves++;
                append_child(children, &count, id);
                if (id == 1) goto generated; /* one P child is a complete N proof */
            }
        }
    }
    for (int k = 0; k < e->k; k++) {
        uint32_t id = proof_key(e->parts, e->count, k);
        append_child(children, &count, id);
        if (id == 1) break;
    }
generated:
    if (edge_used+(size_t)count > edge_capacity) proof_budget();
    e->children = proof_edges+edge_used;
    memcpy(e->children, children, count*sizeof *children);
    e->child_count = (unsigned short)count;
    edge_used += count;
    e->expanded = 1;
    proof_expansions++;
}

static void update(struct proof_node *e) {
    if (e->count == 1) {
        struct query_entry *q = query_slot(e->parts[0]);
        if (q->exact >= 0) {
            e->proof = q->exact == e->k ? INF : 0;
            e->disproof = q->exact == e->k ? 0 : INF;
            return;
        }
        if ((q->excluded >> e->k) & 1) { e->proof = 0; e->disproof = INF; return; }
    }
    e->proof = INF; e->disproof = 0;
    for (int i = 0; i < e->child_count; i++) {
        struct proof_node *child = &proof_table[e->children[i]];
        if (child->disproof < e->proof) e->proof = child->disproof;
        e->disproof = e->disproof >= INF-child->proof ? INF : e->disproof+child->proof;
    }
    if (e->proof == 0) e->disproof = INF;
    if (e->disproof == 0) e->proof = INF;
    if (e->count == 1) {
        struct query_entry *q = query_slot(e->parts[0]);
        if (e->proof == 0) q->excluded |= (Mask)1 << e->k;
        else if (e->disproof == 0) set_exact(q, e->k);
    }
}

static void dfpn(uint32_t id, uint32_t proof_limit, uint32_t disproof_limit) {
    struct proof_node *e = &proof_table[id];
    proof_visits++;
    if (!e->expanded) expand(e);
    update(e);
    while (e->proof < proof_limit && e->disproof < disproof_limit) {
        int selected = -1;
        uint32_t minimum = INF, second = INF;
        for (int i = 0; i < e->child_count; i++) {
            uint32_t d = proof_table[e->children[i]].disproof;
            if (d < minimum) { second = minimum; minimum = d; selected = i; }
            else if (d < second) second = d;
        }
        if (selected < 0) { fprintf(stderr, "dfpn: invalid unsolved node\n"); exit(1); }
        uint32_t child = e->children[selected];
        uint32_t pmax = disproof_limit - e->disproof + proof_table[child].proof;
        uint32_t dmax = second >= proof_limit ? proof_limit : second+1;
        dfpn(child, pmax, dmax);
        update(e);
    }
}

static int proof_sum(Mask *parts, int count, int k) {
    uint32_t id = proof_key(parts, count, k);
    if (id < 2) return id == 1;
    dfpn(id, INF, INF);
    struct proof_node *e = &proof_table[id];
    if (e->proof && e->disproof) { fprintf(stderr, "dfpn: saturated proof numbers; incomplete\n"); exit(3); }
    return e->disproof == 0;
}

static int proof_single(Mask m, int k) { return proof_sum(&m, 1, k); }

int main(int argc, char **argv) {
    const char *env = getenv("STRIP3_DFPN_THRESHOLD");
    if (env) small_threshold = parse_number(env, 0, 20);
    env = getenv("STRIP3_DFPN_BITS");
    int bits = env ? parse_number(env, 8, 22) : 18;
    if (bits < 0 || small_threshold < 0) return 2;
    proof_capacity = (size_t)1 << bits;
    key_capacity = ((size_t)32 << 20)/sizeof *proof_keys;
    edge_capacity = ((size_t)64 << 20)/sizeof *proof_edges;
    proof_table = calloc(proof_capacity, sizeof *proof_table);
    proof_keys = malloc(key_capacity*sizeof *proof_keys);
    proof_edges = malloc(edge_capacity*sizeof *proof_edges);
    if (!proof_table || !proof_keys || !proof_edges) return 3;
    proof_table[0].proof = 0; proof_table[0].disproof = INF;
    proof_table[1].proof = INF; proof_table[1].disproof = 0;
    component_equality_hook = proof_single;
    sum_equality_hook = proof_sum;
    int status = target_cli_main(argc, argv);
    fprintf(stderr, "dfpn proof_memo=%zu expansions=%llu visits=%llu generated_moves=%llu edges=%zu key_bytes=%zu\n",
            proof_entries, (unsigned long long)proof_expansions, (unsigned long long)proof_visits,
            (unsigned long long)proof_generated_moves, edge_used, key_used*sizeof *proof_keys);
    free(proof_table); free(proof_keys); free(proof_edges);
    return status;
}
