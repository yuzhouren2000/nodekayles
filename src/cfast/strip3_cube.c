/* Proof-carrying target queries over independent histories, with no BDD.
 * A cube fixes selected (yes) and unselected (no) history vertices. Facts
 * hold for ALL independent histories in the cube, not just the sample.
 * N: a legal move into a P cube, or membership in a smaller-target P cube.
 * P: every graph move is blocked or enters an N cube; all lower targets N.
 * Facts depend only on earlier facts, so the exported inference DAG has a
 * local checker independent of any game search or symbolic set operations.
 * Five 64-bit words support 100 columns geometrically. This is experimental;
 * the presence of that input size says nothing about solving the full board.
 * Rigidity of equality cubes prevents cross-history generalisation in these
 * rules (notes/37). Retained as a checked counterexperiment, not a new engine.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { WORDS=5, VMAX=300, TARGETS=302 };
typedef struct { uint64_t w[WORDS]; } Bits;
typedef struct { Bits yes, no; } Cube;
typedef struct {
    Cube cube;
    uint32_t next, start, child;
    uint16_t k, move;
    uint8_t truth, kind;
} Fact;
typedef struct { uint32_t child[3], fact; } Trie;
typedef struct { Bits history; uint32_t fact; uint16_t k; } Memo;
static Bits closed[VMAX], adjacent[VMAX];
static Fact *facts;
static Trie *trie;
static Memo *memo;
static uint32_t *edges;
static uint32_t roots[TARGETS][2], heads[TARGETS][2];
static size_t fact_count, fact_limit, trie_count=1, trie_limit, edge_count, edge_limit;
static size_t memo_capacity, memo_count;
static int vertices, columns, lookup_mode=1;
static uint64_t queries, moves, exact_hits, region_hits, matches, trie_visits, literal_sum;
static uint64_t query_limit, next_trace=1024;
static int trace;
static double started;

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
static int has(Bits b,int v) { return (b.w[v/64]>>(v%64))&1; }
static void put(Bits *b,int v) { b->w[v/64]|=UINT64_C(1)<<(v%64); }
static void clear(Bits *b,int v) { b->w[v/64]&=~(UINT64_C(1)<<(v%64)); }
static Bits unite(Bits a,Bits b) { for(int w=0;w<WORDS;w++) a.w[w]|=b.w[w]; return a; }
static int intersects(Bits a,Bits b) { for(int w=0;w<WORDS;w++) if(a.w[w]&b.w[w]) return 1; return 0; }
static int equal_bits(Bits a,Bits b) { for(int w=0;w<WORDS;w++) if(a.w[w]!=b.w[w]) return 0; return 1; }
static int bit_count(Bits b) { int n=0; for(int w=0;w<WORDS;w++) n+=__builtin_popcountll(b.w[w]); return n; }
static int member(Cube c,Bits h) {
    for(int w=0;w<WORDS;w++) if((h.w[w]&c.yes.w[w])!=c.yes.w[w] || (h.w[w]&c.no.w[w])) return 0;
    return 1;
}
static void stats(void) {
    fprintf(stderr,"cube time=%.6f queries=%llu moves=%llu facts=%zu p_facts=",now()-started,
            (unsigned long long)queries,(unsigned long long)moves,fact_count);
    size_t p=0; for(size_t i=1;i<=fact_count;i++) p+=facts[i].truth;
    fprintf(stderr,"%zu edges=%zu literals=%llu exact_hits=%llu region_hits=%llu cube_tests=%llu trie_visits=%llu trie_nodes=%zu memo=%zu lookup=%d\n",
            p,edge_count,(unsigned long long)literal_sum,(unsigned long long)exact_hits,
            (unsigned long long)region_hits,(unsigned long long)matches,
            (unsigned long long)trie_visits,trie_count,memo_count,lookup_mode);
}
static void fail(const char *what) { fprintf(stderr,"cube: %s; incomplete\n",what); exit(3); }
static int independent(Bits h) { for(int v=0;v<vertices;v++) if(has(h,v)&&intersects(h,adjacent[v])) return 0; return 1; }
static Cube normalize(Cube c) {
    if(intersects(c.yes,c.no)||!independent(c.yes)) fail("inconsistent learned condition");
    for(int v=0;v<vertices;v++) if(has(c.yes,v))
        for(int w=0;w<WORDS;w++) c.no.w[w]&=~adjacent[v].w[w];
    return c;
}
static Cube conjunction(Cube a,Cube b) { a.yes=unite(a.yes,b.yes); a.no=unite(a.no,b.no); return a; }
/* Under a legal move at v, v becomes 1 and its neighbours must be 0.
 * These literals need not be fixed by the parent condition. */
static Cube pullback(Cube c,int v) {
    if(has(c.no,v)) fail("false successor explanation");
    clear(&c.yes,v);
    for(int w=0;w<WORDS;w++) c.no.w[w]&=~adjacent[v].w[w];
    return c;
}
static uint64_t hash_history(Bits h,int k) {
    uint64_t x=(uint64_t)k+UINT64_C(0x9e3779b97f4a7c15);
    for(int w=0;w<WORDS;w++) { x^=h.w[w]; x^=x>>30; x*=UINT64_C(0xbf58476d1ce4e5b9); x^=x>>27; x*=UINT64_C(0x94d049bb133111eb); x^=x>>31; }
    return x;
}
static size_t memo_slot(Bits h,int k) {
    size_t i=hash_history(h,k)&(memo_capacity-1);
    while(memo[i].fact && (memo[i].k!=k || !equal_bits(memo[i].history,h))) i=(i+1)&(memo_capacity-1);
    return i;
}
static uint32_t remembered(Bits h,int k) { return memo[memo_slot(h,k)].fact; }
static void remember(Bits h,int k,uint32_t id) {
    size_t i=memo_slot(h,k);
    if(!memo[i].fact) { if(memo_count>=memo_capacity*7/10) fail("exact memo budget"); memo_count++; }
    memo[i]=(Memo){h,id,(uint16_t)k};
}
static uint32_t new_trie(void) { if(trie_count>=trie_limit) fail("condition index budget"); return (uint32_t)trie_count++; }
static void index_fact(uint32_t id) {
    Fact *f=&facts[id];
    f->next=heads[f->k][f->truth]; heads[f->k][f->truth]=id;
    if(lookup_mode!=1) return;
    uint32_t *root=&roots[f->k][f->truth]; if(!*root) *root=new_trie();
    uint32_t at=*root;
    int last=vertices-1;
    while(last>=0&&!has(f->cube.yes,last)&&!has(f->cube.no,last)) last--;
    for(int v=0;v<=last;v++) {
        int b=has(f->cube.yes,v)?1:has(f->cube.no,v)?0:2;
        uint32_t child=trie[at].child[b];
        if(!child) { child=new_trie(); trie[at].child[b]=child; }
        at=child;
    }
    if(!trie[at].fact) trie[at].fact=id;
}
static uint32_t trie_find(uint32_t at,int v,Bits h) {
    if(!at) return 0;
    trie_visits++;
    if(trie[at].fact) return trie[at].fact;
    uint32_t found=trie_find(trie[at].child[2],v+1,h);
    if(found) return found;
    if(v>=vertices) return 0;
    return trie_find(trie[at].child[has(h,v)],v+1,h);
}
static uint32_t region(Bits h,int k,int truth) {
    if(lookup_mode==1) return trie_find(roots[k][truth],0,h);
    if(lookup_mode==2) for(uint32_t i=heads[k][truth];i;i=facts[i].next) { matches++; if(member(facts[i].cube,h)) return i; }
    return 0;
}
static uint32_t learn(Bits h,int k,int truth,int kind,Cube c,int move,uint32_t child,int32_t *refs) {
    c=normalize(c);
    if(!member(c,h)) fail("condition does not cover completed query");
    if(fact_count+1>=fact_limit) fail("proof fact budget");
    uint32_t id=(uint32_t)++fact_count;
    Fact *f=&facts[id];
    *f=(Fact){.cube=c,.k=(uint16_t)k,.truth=(uint8_t)truth,.kind=(uint8_t)kind,.move=(uint16_t)move,.child=child};
    if(truth) {
        size_t count=(size_t)vertices+k;
        if(edge_count+count>edge_limit) fail("proof edge budget");
        f->start=(uint32_t)edge_count;
        for(size_t i=0;i<count;i++) edges[edge_count++]=(uint32_t)refs[i];
    }
    literal_sum+=(unsigned)bit_count(c.yes)+(unsigned)bit_count(c.no);
    index_fact(id); remember(h,k,id);
    return id;
}
static uint32_t query(Bits h,int k) {
    queries++;
    if(trace&&queries>=next_trace) { stats(); next_trace*=2; }
    if(query_limit&&queries>query_limit) fail("query work budget");
    uint32_t id=remembered(h,k);
    if(id) { exact_hits++; return id; }
    for(int truth=1;truth>=0;truth--) if((id=region(h,k,truth))) { region_hits++; remember(h,k,id); return id; }
    for(int j=0;j<k;j++) if((id=region(h,j,1)))
        return learn(h,k,0,'L',facts[id].cube,0,id,NULL);
    Cube cube={0};
    int32_t refs[VMAX*2+1];
    for(int v=0;v<vertices;v++) refs[v]=0;
    /* Middle row first, from the central column towards both ends. */
    for(int pass=0;pass<2;pass++) for(int step=0;step<columns;step++) {
        int c=columns/2+(step%2?-(step+1)/2:step/2);
        for(int r=0;r<3;r++) {
            int v=3*c+r;
            if((r==1)!=(pass==0)) continue;
            if(intersects(h,closed[v])) {
                int u=0; while(!(has(h,u)&&has(closed[v],u))) u++;
                refs[v]=-u-1; put(&cube.yes,u); continue;
            }
            moves++;
            Bits next=h; put(&next,v);
            uint32_t child=query(next,k);
            if(facts[child].truth) {
                Cube explanation=pullback(facts[child].cube,v);
                explanation.no=unite(explanation.no,closed[v]);
                return learn(h,k,0,'M',explanation,v,child,NULL);
            }
            refs[v]=(int32_t)child;
            cube=conjunction(cube,pullback(facts[child].cube,v));
        }
    }
    for(int j=0;j<k;j++) {
        uint32_t child=query(h,j);
        if(facts[child].truth) return learn(h,k,0,'L',facts[child].cube,0,child,NULL);
        refs[vertices+j]=(int32_t)child;
        cube=conjunction(cube,facts[child].cube);
    }
    return learn(h,k,1,'P',cube,0,0,refs);
}
static uint32_t value(Bits h) {
    for(int k=0;k<=vertices;k++) { uint32_t id=query(h,k); if(facts[id].truth) return id; }
    fail("mex bound violated"); return 0;
}
static void write_bits(FILE *f,Bits b) {
    int w=WORDS-1; while(w>0&&!b.w[w]) w--;
    fprintf(f,"%llx",(unsigned long long)b.w[w]);
    while(w--) fprintf(f,"%016llx",(unsigned long long)b.w[w]);
}
static int export_proof(const char *path,uint32_t root) {
    size_t length=strlen(path)+32; char *tmp=malloc(length);
    if(!tmp) return 0;
    snprintf(tmp,length,"%s.tmp",path);
    FILE *f=fopen(tmp,"w"); if(!f) { free(tmp); return 0; }
    fprintf(f,"node-kayles-history-cubes-v1 %d %zu %u\n",columns,fact_count,root);
    for(uint32_t i=1;i<=fact_count;i++) {
        Fact *a=&facts[i]; fprintf(f,"%u %u %u %c ",i,a->k,a->truth,a->kind);
        write_bits(f,a->cube.yes); fputc(' ',f); write_bits(f,a->cube.no);
        if(a->truth) for(int j=0;j<vertices+a->k;j++) fprintf(f," %d",(int32_t)edges[a->start+j]);
        else if(a->kind=='M') fprintf(f," %u %u",a->move,a->child);
        else fprintf(f," %u",a->child);
        fputc('\n',f);
    }
    int ok=!ferror(f); if(fclose(f)) ok=0;
    if(ok&&rename(tmp,path)) ok=0;
    if(!ok) remove(tmp);
    free(tmp); return ok;
}
static int number(const char *s,int lo,int hi) {
    char *end; long value=strtol(s,&end,10); return !*s||*end||value<lo||value>hi?-1:(int)value;
}
static int parse_bits(const char *s,Bits *h) {
    memset(h,0,sizeof *h);
    if(!*s||strlen(s)>(size_t)(vertices+3)/4) return 0;
    for(const char *p=s;*p;p++) {
        int d=*p>='0'&&*p<='9'?*p-'0':*p>='a'&&*p<='f'?*p-'a'+10:*p>='A'&&*p<='F'?*p-'A'+10:-1;
        if(d<0) return 0;
        for(int w=WORDS-1;w>=0;w--) h->w[w]=(h->w[w]<<4)|(w?h->w[w-1]>>60:(unsigned)d);
    }
    for(int v=vertices;v<WORDS*64;v++) if(has(*h,v)) return 0;
    return independent(*h);
}
int main(int argc,char **argv) {
    const char *certificate=NULL;
    if(argc==5&&!strcmp(argv[3],"--certificate")) { certificate=argv[4]; argc=3; }
    if(argc!=3) return 2;
    int board=!strcmp(argv[1],"board"), histories=!strcmp(argv[1],"histories"), qs=!strcmp(argv[1],"queries");
    columns=number(argv[2],1,100); if(columns<0||(!board&&!histories&&!qs)||(certificate&&!board)) return 2;
    vertices=3*columns;
    const char *s=getenv("STRIP3_CUBE_LOOKUP"); if(s) lookup_mode=number(s,0,2); if(lookup_mode<0) return 2;
    s=getenv("STRIP3_CUBE_BITS"); int bits=s?number(s,8,20):18; if(bits<0) return 2;
    s=getenv("STRIP3_CUBE_TRACE"); if(s) trace=number(s,0,1); if(trace<0) return 2;
    s=getenv("STRIP3_CUBE_QUERY_LIMIT"); if(s) { char *end; if(*s=='-') return 2; query_limit=strtoull(s,&end,10); if(!*s||*end||!query_limit) return 2; }
    fact_limit=(size_t)1<<bits; memo_capacity=fact_limit*2; trie_limit=fact_limit*8; edge_limit=fact_limit*16;
    facts=calloc(fact_limit,sizeof *facts); memo=calloc(memo_capacity,sizeof *memo);
    trie=calloc(trie_limit,sizeof *trie); edges=calloc(edge_limit,sizeof *edges);
    if(!facts||!memo||!trie||!edges) { fprintf(stderr,"cube: allocation failed; incomplete\n"); return 3; }
    for(int v=0;v<vertices;v++) for(int u=0;u<vertices;u++) {
        int distance=abs(v/3-u/3)+abs(v%3-u%3);
        if(distance<=1) put(&closed[v],u);
        if(distance==1) put(&adjacent[v],u);
    }
    started=now(); atexit(stats);
    if(board) {
        Bits h={0}; uint32_t root=value(h);
        if(certificate&&!export_proof(certificate,root)) fail("proof export failed");
        printf("3\t%d\t%u\n",columns,facts[root].k);
    } else {
        char line[1024],hex[100],extra; int k;
        while(fgets(line,sizeof line,stdin)) {
            Bits h;
            if(qs) { if(sscanf(line,"%99s %d %c",hex,&k,&extra)!=2||k<0||k>vertices||!parse_bits(hex,&h)) return 2; }
            else if(sscanf(line,"%99s %c",hex,&extra)!=1||!parse_bits(hex,&h)) return 2;
            uint32_t root=qs?query(h,k):value(h);
            printf("%u\n",qs?facts[root].truth:facts[root].k); fflush(stdout);
        }
        if(ferror(stdin)) return 2;
    }
    return 0;
}
