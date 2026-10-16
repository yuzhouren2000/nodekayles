/* Export proofs from the fast target search's existing knowledge only.
 * No SG oracle or new equality search runs during certificate extraction.
 * CLI: board N --certificate PATH. All inherited search controls still apply.
 */
#ifdef STRIP3_TARGET_PROOF_LAZY
#define STRIP3_LAZY_EXTERNAL_MAIN
#include "strip3_lazy.c"
#else
#define TARGET_ENTRY_POINT target_cli_main
#include "strip3_target.c"
#endif
#define STRIP3_PROOF_EXTERNAL_BASE
#define STRIP3_PROOF_EXTERNAL_MAIN
#include "strip3_proof_cost.c"

static const char *certificate_path;
static double overall_start, search_time, export_time;
static Mask exported_cost;
static int exported_n, exported_g;
static size_t search_entries;

static struct query_entry *known_entry(Mask m) {
    size_t i = hash(m)&(query_capacity-1);
    while (query_table[i].key && query_table[i].key != m) i = (i+1)&(query_capacity-1);
    return query_table[i].key ? &query_table[i] : NULL;
}

static int known_value(Mask m) {
    struct query_entry *entry = known_entry(m);
    if (entry && entry->exact >= 0) return entry->exact;
    int rule = axiom(m), nodes = order(m);
    if (rule == R_ROTATION) return 0;
    if (rule == R_PATH) return path_g[nodes];
    if (rule == R_CYCLE) return path_g[nodes-3] == 0;
    return -1;
}

static int known_truth(Mask m, int k) {
    if (k > order(m)) return 0;
    int g = known_value(m);
    if (g >= 0) return g == k;
    struct query_entry *entry = known_entry(m);
    return entry && (entry->excluded >> k & 1) ? 0 : -1;
}

static void export_completed(Mask root, int g) {
    double now = seconds(); search_time = now-overall_start;
    search_entries = query_entries;
    uint64_t before_queries = queries, before_moves = query_moves;
    proof_value_hook = known_value; proof_truth_hook = known_truth;
    exported_cost = tree_cost(root, g);
    if (exported_cost == COST_INF) {
        fprintf(stderr, "target-proof: missing or oversized proof in available knowledge; incomplete\n"); exit(3);
    }
    emit_certificate(certificate_path, width(root), root, g, exported_cost);
    if (tst.calls || tst.moves || entries || query_entries != search_entries ||
            queries != before_queries || query_moves != before_moves) {
        fprintf(stderr, "target-proof: extraction performed new game search; incomplete\n"); exit(4);
    }
    export_time = seconds()-now; exported_n = width(root); exported_g = g;
}

int main(int argc, char **argv) {
    if (argc != 5 || strcmp(argv[1], "board") || strcmp(argv[3], "--certificate")) {
        fprintf(stderr, "usage: strip3_target_proof board N --certificate PATH\n"); return 2;
    }
    const char *env = getenv("STRIP3_PROOF_BITS"); int bits = env ? integer(env, 8, 24) : 20;
    if (bits < 0) return 2;
    env = getenv("STRIP3_PROOF_FIRST");
    if (env && (first_policy = integer(env, 0, 1)) < 0) return 2;
    proof_cap = (size_t)1 << bits; proofs = calloc(proof_cap, sizeof *proofs);
    if (!proofs) { fprintf(stderr, "target-proof: allocation failed; incomplete\n"); return 3; }
    certificate_path = argv[4]; completed_board_hook = export_completed;
#ifdef STRIP3_TARGET_PROOF_LAZY
    int configured = configure_lazy();
    if (configured) return configured;
#endif
    overall_start = seconds();
    int status = target_cli_main(3, argv);
    if (status) return status;
    fprintf(stdout, "{\"n\":%d,\"SG\":%d,\"policy\":\"%s\",\"tree_cost\":\"", exported_n, exported_g, policy_name());
    decimal(stdout, exported_cost);
    fprintf(stdout, "\",\"search_states\":%zu,\"proof_states\":%zu,\"certificate_queries\":%zu,\"sum_rules\":%llu,\"edges\":%llu,\"considered_moves\":%llu,\"sum_evaluations\":%llu,\"oracle_calls\":%llu,\"search_s\":%.6f,\"export_s\":%.6f,\"total_s\":%.6f}\n",
            search_entries, proof_count, selected_count, (unsigned long long)emitted_sum_rules,
            (unsigned long long)emitted_edges, (unsigned long long)considered_moves,
            (unsigned long long)sum_evaluations, (unsigned long long)tst.calls,
            search_time, export_time, seconds()-overall_start);
    return 0;
}
