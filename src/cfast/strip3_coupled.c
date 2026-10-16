/* Experimental P/N search on a disjoint sum of unresolved components and *k.
 * This removes the requirement to evaluate every smaller component exactly.
 * It is exact, but the tuple memo can grow faster than a component SG memo.
 * Reuse geometry and the tested single-component query engine; a sum query
 * explores graph moves and Nim moves directly, with an exact tuple key.
 */
#define TARGET_ENTRY_POINT target_cli_main
#include "strip3_target.c"

struct sum_entry { Mask *parts; uint64_t hash; unsigned short count; unsigned char k, outcome; };
static struct sum_entry *sum_table;
static Mask *sum_keys;
static size_t sum_capacity, sum_entries, sum_key_capacity, sum_key_used;
static uint64_t sum_queries, sum_moves;
static int exact_threshold = 12, nim_first;

static int coupled_equality(Mask *parts, int count, int k);

static void sum_budget(void) {
    fprintf(stderr, "coupled: memo budget reached; query incomplete; sum_memo=%zu sum_moves=%llu\n",
            sum_entries, (unsigned long long)sum_moves);
    exit(3);
}

static int known_component(Mask m) {
    size_t i = hash(m) & (query_capacity - 1);
    while (query_table[i].key) {
        if (query_table[i].key == m && query_table[i].exact >= 0) return query_table[i].exact;
        if (query_table[i].key == m) break;
        i = (i + 1) & (query_capacity - 1);
    }
    int n = width(m), h = (m & rowmask[2]) ? 3 : ((m & rowmask[1]) ? 2 : 1);
    Mask both = reverse_bits(m, 3*n) >> (3-h);
    if (m == both) {
        int safe = 1;
        for (Mask x = m; x && safe; x &= x-1) {
            int v = lowbit(x), partner = 3*(n-1-v/3)+h-1-v%3;
            if (partner == v || (neighbours((Mask)1 << v) & ((Mask)1 << partner))) safe = 0;
        }
        if (safe) return 0;
    }
    int degree_max = 0, endpoints = 0;
    for (Mask x = m; x; x &= x-1) {
        int degree = order(neighbours((Mask)1 << lowbit(x)) & m);
        if (degree > degree_max) degree_max = degree;
        if (degree < 2) endpoints++;
        if (degree_max > 2) break;
    }
    if (degree_max <= 2)
        return endpoints || order(m) == 1 ? path_g[order(m)] : (path_g[order(m)-3] == 0);
    if (order(m) <= exact_threshold) return exact_value(m);
    return -1;
}

static void add_part(Mask *parts, int *count, Mask key) {
    for (int i = 0; i < *count; i++) if (parts[i] == key) {
        parts[i] = parts[--*count];
        return;
    }
    parts[(*count)++] = key;
}

static uint64_t tuple_hash(Mask *parts, int count, int k) {
    uint64_t h = (uint64_t)(count*128+k) * 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < count; i++) h = (h ^ hash(parts[i])) * 0xbf58476d1ce4e5b9ULL;
    return h ^ (h >> 31);
}

static struct sum_entry *sum_slot(Mask *parts, int count, int k) {
    uint64_t h = tuple_hash(parts, count, k);
    size_t i = h & (sum_capacity-1);
    while (sum_table[i].count) {
        struct sum_entry *e = &sum_table[i];
        if (e->hash == h && e->count == count && e->k == k && !memcmp(e->parts, parts, count*sizeof *parts))
            return e;  /* hashes alone never establish key equality */
        i = (i+1) & (sum_capacity-1);
    }
    if (sum_entries >= sum_capacity*7/10 || sum_key_used+(size_t)count > sum_key_capacity) sum_budget();
    struct sum_entry *e = &sum_table[i];
    e->parts = sum_keys+sum_key_used;
    memcpy(e->parts, parts, count*sizeof *parts);
    sum_key_used += count;
    e->hash = h; e->count = (unsigned short)count; e->k = (unsigned char)k;
    sum_entries++;
    return e;
}

static int coupled_equality(Mask *input, int count, int k) {
    sum_queries++;
    Mask parts[MAXCELLS];
    int kept = 0, total_nodes = 0;
    for (int i = 0; i < count; i++) {
        int g = known_component(input[i]);
        if (g >= 0) k ^= g;
        else { parts[kept++] = input[i]; total_nodes += order(input[i]); }
    }
    count = kept;
    if (!count) return k == 0;
    if (k > total_nodes) return 0; /* xor of nimbers <= sum of component orders */
    if (count == 1) return equality_component(parts[0], k);
    /* Stable tuple key independent of component order. */
    for (int i = 1; i < count; i++) {
        Mask x = parts[i]; int j = i;
        while (j && parts[j-1] > x) { parts[j] = parts[j-1]; j--; }
        parts[j] = x;
    }
    struct sum_entry *entry = sum_slot(parts, count, k);
    if (entry->outcome) return entry->outcome == 2;
    if (nim_first) for (int j = 0; j < k; j++) if (coupled_equality(parts, count, j)) {
        entry->outcome = 1;
        return 0;
    }
    int indexes[MAXCELLS];
    for (int i = 0; i < count; i++) indexes[i] = i;
    for (int i = 1; i < count; i++) {
        int x = indexes[i], j = i;
        while (j && order(parts[indexes[j-1]]) < order(parts[x])) { indexes[j] = indexes[j-1]; j--; }
        indexes[j] = x;
    }
    for (int a = 0; a < count; a++) {
        int selected = indexes[a];
        Mask m = parts[selected];
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
                Mask next[MAXCELLS], children[MAXCELLS]; int nc = 0;
                for (int i = 0; i < count; i++) if (i != selected) add_part(next, &nc, parts[i]);
                Mask bit = (Mask)1 << v;
                int cn = split_position(m & ~(bit | neighbours(bit)), children);
                for (int i = 0; i < cn; i++) add_part(next, &nc, children[i]);
                sum_moves++;
                if (coupled_equality(next, nc, k)) { entry->outcome = 1; return 0; }
            }
        }
    }
    if (!nim_first) for (int j = 0; j < k; j++) if (coupled_equality(parts, count, j)) {
        entry->outcome = 1;
        return 0;
    }
    entry->outcome = 2;
    return 1;
}

int main(int argc, char **argv) {
    const char *env = getenv("STRIP3_COUPLED_THRESHOLD");
    if (env) exact_threshold = parse_number(env, 0, 126);
    env = getenv("STRIP3_COUPLED_NIM_FIRST");
    if (env) nim_first = parse_number(env, 0, 1);
    env = getenv("STRIP3_COUPLED_BITS");
    int bits = env ? parse_number(env, 8, 22) : 18;
    if (exact_threshold < 0 || nim_first < 0 || bits < 0) return 2;
    sum_capacity = (size_t)1 << bits;
    sum_key_capacity = ((size_t)64 << 20) / sizeof *sum_keys;
    sum_table = calloc(sum_capacity, sizeof *sum_table);
    sum_keys = malloc(sum_key_capacity*sizeof *sum_keys);
    if (!sum_table || !sum_keys) { fprintf(stderr, "coupled: allocation failed\n"); return 3; }
    sum_equality_hook = coupled_equality;
    int result = target_cli_main(argc, argv);
    fprintf(stderr, "coupled sum_memo=%zu sum_queries=%llu sum_moves=%llu key_bytes=%zu threshold=%d nim_first=%d\n",
            sum_entries, (unsigned long long)sum_queries, (unsigned long long)sum_moves,
            sum_key_used*sizeof *sum_keys, exact_threshold, nim_first);
    free(sum_table); free(sum_keys);
    return result;
}
