/* Bounded experiment: replace arbitrary caterpillar components by an exact
 * O(p^3)-time, O(p^2)-space interval DP on their p spine vertices.
 * Recognition checks the graph itself, including for unreachable input masks.
 * This is not an algorithm for thick full strips or for 3x100.
 */
#define TARGET_ENTRY_POINT target_cli_main
#include "strip3_target.c"

static uint64_t cat_checks, cat_squares, cat_scanned, cat_hits;
static uint64_t cat_intervals, cat_options;
static int cat_max_spine;

static int caterpillar_value(Mask m) {
    cat_checks++;
    /* A tree contains no unit square; reject thick masks without a scan. */
    if (m & (m >> 1) & (m >> 3) & (m >> 4) & (rowmask[0] | rowmask[1])) {
        cat_squares++;
        return -1;
    }
    int nodes = order(m), degree_sum = 0;
    Mask spine = 0;
    for (Mask x = m; x; x &= x-1) {
        int v = lowbit(x), d = order(neighbours((Mask)1 << v) & m);
        cat_scanned++;
        degree_sum += d;
        if (d >= 2) spine |= (Mask)1 << v;
    }
    if (!nodes || degree_sum != 2*(nodes-1)) return -1;
    /* The hook receives connected components. Check connectedness too so the
     * evaluator is safe when tested directly on arbitrary input masks. */
    Mask reached = (Mask)1 << lowbit(m), frontier = reached;
    while (frontier) { frontier = neighbours(frontier) & m & ~reached; reached |= frontier; }
    if (reached != m) return -1;
    if (!spine) { cat_hits++; return 1; }  /* K1 or K2 */
    int p = order(spine), start = -1;
    for (Mask x = spine; x; x &= x-1) {
        int v = lowbit(x), d = order(neighbours((Mask)1 << v) & spine);
        cat_scanned++;
        if (d > 2) return -1;
        if (d < 2) start = v;
    }
    if (start < 0) return -1;
    unsigned char leaves[MAXCELLS], dp[MAXCELLS][MAXCELLS];
    Mask remaining = spine;
    int v = start;
    for (int i = 0; i < p; i++) {
        Mask bit = (Mask)1 << v;
        leaves[i] = (unsigned char)order(neighbours(bit) & m & ~spine);
        remaining &= ~bit;
        Mask next = neighbours(bit) & remaining;
        if (i+1 < p) {
            if (!next || (next & (next-1))) return -1;
            v = lowbit(next);
        } else if (remaining) return -1;
    }
    /* Every interval keeps its original attached leaves. Playing a spine
     * vertex frees the leaves on the two deleted neighbouring spine vertices;
     * playing a leaf frees its siblings. Isolated vertices xor by parity. */
    for (int length = 1; length <= p; length++) for (int a = 0; a+length <= p; a++) {
        int b = a+length-1;
        Mask seen = 0;
        cat_intervals++;
        for (int i = a; i <= b; i++) {
            int g = (i >= a+2 ? dp[a][i-2] : 0) ^ (i+2 <= b ? dp[i+2][b] : 0);
            if (i > a) g ^= leaves[i-1] & 1;
            if (i < b) g ^= leaves[i+1] & 1;
            seen |= (Mask)1 << g;
            cat_options++;
            if (leaves[i]) {
                g = (i > a ? dp[a][i-1] : 0) ^ (i < b ? dp[i+1][b] : 0)
                    ^ ((leaves[i]-1) & 1);
                seen |= (Mask)1 << g;
                cat_options++;
            }
        }
        int g = 0;
        while ((seen >> g) & 1) g++;
        dp[a][b] = (unsigned char)g;
    }
    cat_hits++;
    if (p > cat_max_spine) cat_max_spine = p;
    return dp[0][p-1];
}

#ifndef STRIP3_CATERPILLAR_EXTERNAL_MAIN
int main(int argc, char **argv) {
    const char *env = getenv("STRIP3_CATERPILLAR");
    int enabled = env ? parse_number(env, 0, 1) : 1;
    if (enabled < 0) { fprintf(stderr, "STRIP3_CATERPILLAR must be 0 or 1\n"); return 2; }
    if (enabled) component_value_hook = caterpillar_value;
    int status = target_cli_main(argc, argv);
    fprintf(stderr, "caterpillar enabled=%d cat_checks=%llu cat_squares=%llu cat_scanned=%llu cat_hits=%llu cat_max_spine=%d cat_intervals=%llu cat_options=%llu\n",
            enabled, (unsigned long long)cat_checks, (unsigned long long)cat_squares,
            (unsigned long long)cat_scanned, (unsigned long long)cat_hits, cat_max_spine,
            (unsigned long long)cat_intervals, (unsigned long long)cat_options);
    return status;
}
#endif
