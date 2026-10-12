/* kayles.c - fast exact Sprague-Grundy solver for Node Kayles on grid boards.
 *
 * Same algorithm as the Python reference in src/nodekayles/solver.py:
 *   - a position is the set of playable cells; the value is the xor of the
 *     values of its connected components (Sprague-Grundy decomposition)
 *   - moves: place a stone at a playable cell v, deleting v and its four grid
 *     neighbours from the playable set
 *   - components are memoised under a canonical form (bounding box, transposed
 *     so that height <= width, minimised over the symmetries of the box)
 *
 * Boards up to MAXBITS = 384 cells are supported (e.g. 3x128, 8x48, 12x12).
 *
 * Usage:
 *   kayles strip M NMAX        print "M n G" for n = 1..NMAX (shares one memo)
 *   kayles board M N           print the Grundy value of the M x N board
 *   kayles table MMAX NMAX     print "m n G" for all 1 <= m <= n <= NMAX
 *
 * Env:
 *   KAYLES_HT_MAXBITS  max log2 of hash capacity (default 25)
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define WORDS 6
#define MAXBITS (64 * WORDS)

typedef struct {
    uint16_t hw; /* h << 8 | w */
    uint64_t m[WORDS];
} Key;

typedef struct {
    Key k;
    int32_t val;
    uint8_t used;
} Entry;

static Entry *ht = NULL;
static size_t ht_cap = 0, ht_cnt = 0, ht_maxbits = 25;

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static uint64_t hashkey(Key k) {
    uint64_t h = 1469598103934665603ULL ^ ((uint64_t)k.hw * 0x9E3779B97F4A7C15ULL);
    for (int i = 0; i < WORDS; i++) {
        h ^= k.m[i];
        h *= 1099511628211ULL;
        h ^= h >> 29;
    }
    return h;
}

static void ht_init(size_t bits) {
    ht_cap = (size_t)1 << bits;
    ht_cnt = 0;
    ht = (Entry *)calloc(ht_cap, sizeof(Entry));
    if (!ht) {
        fprintf(stderr, "kayles: out of memory for hash table (%zu entries)\n", ht_cap);
        exit(1);
    }
}

static void ht_insert_raw(Key k, int32_t val) {
    size_t i = hashkey(k) & (ht_cap - 1);
    while (ht[i].used) i = (i + 1) & (ht_cap - 1);
    ht[i].k = k;
    ht[i].val = val;
    ht[i].used = 1;
}

static void ht_grow(void) {
    if (ht_cap == ((size_t)1 << ht_maxbits)) {
        return; /* keep probing; slower but no OOM */
    }
    size_t newbits = 0;
    for (size_t c = ht_cap; c; c >>= 1) newbits++;
    if (newbits >= ht_maxbits) newbits = ht_maxbits;
    Entry *old = ht;
    size_t oldcap = ht_cap, oldcnt = ht_cnt;
    ht_init(newbits);
    for (size_t i = 0; i < oldcap; i++)
        if (old[i].used) ht_insert_raw(old[i].k, old[i].val);
    free(old);
    ht_cnt = oldcnt;
}

static int ht_get(Key k, int32_t *val) {
    size_t i = hashkey(k) & (ht_cap - 1);
    while (ht[i].used) {
        if (ht[i].k.hw == k.hw && memcmp(ht[i].k.m, k.m, sizeof k.m) == 0) {
            *val = ht[i].val;
            return 1;
        }
        i = (i + 1) & (ht_cap - 1);
    }
    return 0;
}

static void ht_put(Key k, int32_t val) {
    if ((ht_cnt + 1) * 10 >= ht_cap * 7) ht_grow();
    ht_insert_raw(k, val);
    ht_cnt++;
}

/* neighbour bitmask of cell v in an h x w box */
static void nbrs_of(int h, int w, int v, uint64_t *out) {
    memset(out, 0, WORDS * sizeof(uint64_t));
    int r = v / w, c = v % w;
    if (r > 0) out[(v - w) >> 6] |= 1ULL << ((v - w) & 63);
    if (r + 1 < h) out[(v + w) >> 6] |= 1ULL << ((v + w) & 63);
    if (c > 0) out[(v - 1) >> 6] |= 1ULL << ((v - 1) & 63);
    if (c + 1 < w) out[(v + 1) >> 6] |= 1ULL << ((v + 1) & 63);
}

static void transpose_bits(int h, int w, const uint64_t *in, uint64_t *out) {
    memset(out, 0, WORDS * sizeof(uint64_t));
    for (int wi = 0; wi < WORDS; wi++) {
        uint64_t x = in[wi];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int v = wi * 64 + b;
            int r = v / w, c = v % w;
            int nv = c * h + r;
            out[nv >> 6] |= 1ULL << (nv & 63);
        }
    }
}

/* canonical form: bounding box, height <= width, min over box symmetries */
static Key canonical(int h, int w, const uint64_t *mask) {
    int rmin = h, rmax = -1, cmin = w, cmax = -1;
    for (int wi = 0; wi < WORDS; wi++) {
        uint64_t x = mask[wi];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int v = wi * 64 + b;
            int r = v / w, c = v % w;
            if (r < rmin) rmin = r;
            if (r > rmax) rmax = r;
            if (c < cmin) cmin = c;
            if (c > cmax) cmax = c;
        }
    }
    Key out;
    memset(&out, 0, sizeof out);
    if (rmax < 0) return out;
    int bh = rmax - rmin + 1, bw = cmax - cmin + 1;
    uint64_t a[WORDS], cur[WORDS], tmp[WORDS];
    memset(a, 0, sizeof a);
    for (int wi = 0; wi < WORDS; wi++) {
        uint64_t x = mask[wi];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int v = wi * 64 + b;
            int r = v / w, c = v % w;
            int nv = (r - rmin) * bw + (c - cmin);
            a[nv >> 6] |= 1ULL << (nv & 63);
        }
    }
    int H = bh, W = bw;
    memcpy(cur, a, sizeof a);
    if (bh > bw) {
        transpose_bits(bh, bw, a, tmp);
        memcpy(cur, tmp, sizeof tmp);
        H = bw;
        W = bh;
    }
    uint64_t best[WORDS];
    memcpy(best, cur, sizeof cur);
    int ntransforms = (H == W) ? 8 : 4;
    for (int t = 0; t < ntransforms; t++) {
        memset(tmp, 0, sizeof tmp);
        for (int wi = 0; wi < WORDS; wi++) {
            uint64_t x = cur[wi];
            while (x) {
                int b = __builtin_ctzll(x);
                x &= x - 1;
                int v = wi * 64 + b;
                int r = v / W, c = v % W;
                int nr, nc;
                switch (t) {
                    case 0: nr = r; nc = c; break;
                    case 1: nr = H - 1 - r; nc = c; break;
                    case 2: nr = r; nc = W - 1 - c; break;
                    case 3: nr = H - 1 - r; nc = W - 1 - c; break;
                    case 4: nr = c; nc = r; break;
                    case 5: nr = c; nc = H - 1 - r; break;
                    case 6: nr = W - 1 - c; nc = r; break;
                    default: nr = W - 1 - c; nc = H - 1 - r; break;
                }
                int nv = nr * W + nc;
                tmp[nv >> 6] |= 1ULL << (nv & 63);
            }
        }
        if (memcmp(tmp, best, sizeof best) < 0) memcpy(best, tmp, sizeof best);
    }
    out.hw = (uint16_t)((H << 8) | W);
    memcpy(out.m, best, sizeof best);
    return out;
}

static int value_of(Key k);
static int xor_components(int h, int w, const uint64_t *mask);

/* apply transform t of the h x w box to a single cell index */
static int map_cell(int h, int w, int t, int v) {
    int r = v / w, c = v % w, nr, nc;
    switch (t) {
        case 0: nr = r; nc = c; break;
        case 1: nr = h - 1 - r; nc = c; break;
        case 2: nr = r; nc = w - 1 - c; break;
        case 3: nr = h - 1 - r; nc = w - 1 - c; break;
        case 4: nr = c; nc = r; break;
        case 5: nr = c; nc = h - 1 - r; break;
        case 6: nr = w - 1 - c; nc = r; break;
        default: nr = w - 1 - c; nc = h - 1 - r; break;
    }
    return nr * w + nc;
}

static int value_of(Key k) {
    int32_t val;
    if (ht_get(k, &val)) return val;
    int h = k.hw >> 8, w = k.hw & 0xff;
    /* subgroup of box symmetries preserving the position, for move dedup */
    int sub[8], nsub = 1;
    sub[0] = 0;
    int ntrans = (h == w) ? 8 : 4;
    for (int t = 1; t < ntrans; t++) {
        uint64_t img[WORDS];
        memset(img, 0, sizeof img);
        for (int wi = 0; wi < WORDS; wi++) {
            uint64_t x = k.m[wi];
            while (x) {
                int b = __builtin_ctzll(x);
                x &= x - 1;
                int nv = map_cell(h, w, t, wi * 64 + b);
                img[nv >> 6] |= 1ULL << (nv & 63);
            }
        }
        if (memcmp(img, k.m, sizeof img) == 0) sub[nsub++] = t;
    }
    uint64_t seen[2] = {0, 0};
    uint64_t tmp[WORDS], nb[WORDS];
    for (int wi = 0; wi < WORDS; wi++) {
        uint64_t x = k.m[wi];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int v = wi * 64 + b;
            int dedup = 0;
            for (int si = 1; si < nsub && !dedup; si++) {
                int u = map_cell(h, w, sub[si], v);
                if (u < v) dedup = 1;
            }
            if (dedup) continue;
            memcpy(tmp, k.m, sizeof tmp);
            tmp[v >> 6] &= ~(1ULL << (v & 63));
            nbrs_of(h, w, v, nb);
            for (int i = 0; i < WORDS; i++) tmp[i] &= ~nb[i];
            int s = xor_components(h, w, tmp);
            if (s < 128) seen[s >> 6] |= 1ULL << (s & 63);
        }
    }
    int g = 0;
    while (g < 128 && ((seen[g >> 6] >> (g & 63)) & 1)) g++;
    ht_put(k, g);
    return g;
}

static int xor_components(int h, int w, const uint64_t *mask) {
    uint64_t rem[WORDS], vis[WORDS];
    memcpy(rem, mask, sizeof rem);
    int total = 0;
    int stack[MAXBITS];
    for (;;) {
        int v = -1;
        for (int i = 0; i < WORDS; i++)
            if (rem[i]) { v = i * 64 + __builtin_ctzll(rem[i]); break; }
        if (v < 0) break;
        memset(vis, 0, sizeof vis);
        int sp = 0;
        vis[v >> 6] |= 1ULL << (v & 63);
        stack[sp++] = v;
        while (sp) {
            int u = stack[--sp];
            uint64_t nb[WORDS];
            nbrs_of(h, w, u, nb);
            for (int i = 0; i < WORDS; i++) {
                uint64_t add = nb[i] & rem[i] & ~vis[i];
                while (add) {
                    int b = __builtin_ctzll(add);
                    add &= add - 1;
                    int nv = i * 64 + b;
                    vis[nv >> 6] |= 1ULL << (nv & 63);
                    stack[sp++] = nv;
                }
            }
        }
        for (int i = 0; i < WORDS; i++) rem[i] &= ~vis[i];
        Key ck = canonical(h, w, vis);
        total ^= value_of(ck);
    }
    return total;
}

static void build_full(int h, int w, uint64_t *mask) {
    memset(mask, 0, WORDS * sizeof(uint64_t));
    for (int v = 0; v < h * w; v++) mask[v >> 6] |= 1ULL << (v & 63);
}

static int board_value(int h, int w) {
    uint64_t full[WORDS];
    build_full(h, w, full);
    return xor_components(h, w, full);
}

/* Grundy value of the board after a first move at cell v (closed neighbourhood
 * of v deleted). */
static int value_after_move(int h, int w, int v) {
    uint64_t full[WORDS], nb[WORDS];
    build_full(h, w, full);
    nbrs_of(h, w, v, nb);
    full[v >> 6] &= ~(1ULL << (v & 63));
    for (int i = 0; i < WORDS; i++) full[i] &= ~nb[i];
    return xor_components(h, w, full);
}

/* one representative per orbit of cells under the box symmetry group */
static int is_representative(int h, int w, int v) {
    int ntrans = (h == w) ? 8 : 4;
    for (int t = 1; t < ntrans; t++)
        if (map_cell(h, w, t, v) < v) return 0;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: kayles strip M NMAX | board M N | table MMAX NMAX\n");
        return 2;
    }
    const char *env = getenv("KAYLES_HT_MAXBITS");
    if (env) ht_maxbits = (size_t)atoi(env);
    if (ht_maxbits < 16) ht_maxbits = 16;
    double t0 = now_sec();

    if (!strcmp(argv[1], "strip") && argc == 4) {
        int m = atoi(argv[2]), nmax = atoi(argv[3]);
        if (m <= 0 || nmax <= 0 || m * nmax > MAXBITS) {
            fprintf(stderr, "kayles: board too large (max %d cells)\n", MAXBITS);
            return 2;
        }
        ht_init(20);
        for (int n = 1; n <= nmax; n++) {
            double tn = now_sec();
            int g = board_value(m, n);
            printf("%d\t%d\t%d\n", m, n, g);
            fflush(stdout);
            fprintf(stderr, "[strip %dx%d] G=%d  %6.2fs  memo=%zu\n",
                    m, n, g, now_sec() - tn, ht_cnt);
        }
    } else if (!strcmp(argv[1], "board") && argc == 4) {
        int m = atoi(argv[2]), n = atoi(argv[3]);
        if (m * n > MAXBITS) {
            fprintf(stderr, "kayles: board too large (max %d cells)\n", MAXBITS);
            return 2;
        }
        ht_init(20);
        int g = board_value(m, n);
        printf("%d\t%d\t%d\n", m, n, g);
        fflush(stdout);
    } else if (!strcmp(argv[1], "wmoves") && argc == 4) {
        int m = atoi(argv[2]), n = atoi(argv[3]);
        if (m * n > MAXBITS) {
            fprintf(stderr, "kayles: board too large (max %d cells)\n", MAXBITS);
            return 2;
        }
        ht_init(20);
        int g = board_value(m, n);
        printf("# board %dx%d grundy %d\n", m, n, g);
        for (int v = 0; v < m * n; v++) {
            if (!is_representative(m, n, v)) continue;
            int s = value_after_move(m, n, v);
            printf("%d\t%d\t%d\t%s\n", v / n, v % n, s, s == 0 ? "WIN" : "");
            fflush(stdout);
        }
        fprintf(stderr, "wmoves %dx%d done, memo=%zu\n", m, n, ht_cnt);
    } else if (!strcmp(argv[1], "table") && argc == 4) {
        int mmax = atoi(argv[2]), nmax = atoi(argv[3]);
        ht_init(20);
        for (int m = 1; m <= mmax; m++) {
            for (int n = m; n <= nmax; n++) {
                if (m * n > MAXBITS) continue;
                double tn = now_sec();
                int g = board_value(m, n);
                printf("%d\t%d\t%d\n", m, n, g);
                fflush(stdout);
                fprintf(stderr, "[table %dx%d] G=%d  %6.2fs  memo=%zu\n",
                        m, n, g, now_sec() - tn, ht_cnt);
            }
        }
    } else {
        fprintf(stderr, "usage: kayles strip M NMAX | board M N | table MMAX NMAX\n");
        return 2;
    }
    fprintf(stderr, "done in %.1fs, memo entries=%zu, cap=%zu\n",
            now_sec() - t0, ht_cnt, ht_cap);
    return 0;
}
