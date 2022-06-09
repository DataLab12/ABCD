/*
Code for Jelena's plus/minus research project

Copyright (c) 2020 Texas State University. All rights reserved.

Redistribution in source or binary form, with or without modification,
is *not* permitted. Use in source and binary forms, with or without
modification, is only permitted for research use at Texas State University
in Martin Burtscher's, Jelena Tesic's, and Lucas Rusnak's research groups.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Author: Martin Burtscher
*/

#include <cstdio>
#include <climits>
#include <algorithm>
#include <set>
#include <map>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <random>
#include "Printer.h"
#include "DataStructures.h"

// code does not fully adhere to latest OpenMP specification

static const bool verify = true;  // set to false for better performance
static const bool pick_random_roots = false;  // pick tree roots randomly (or from highest to lowest degree)

#define DEBUG 1
#define TEST 1 //for use with test datasets

struct EdgeInfo {
    int beg;  // beginning of range (shifted by 1) | is range inverted or not
    int end;  // end of range (shifted by 1) | plus or minus (1 = minus, 0 = plus or zero)
};

//struct Graph {
//    int nodes;
//    int edges;
//    int* nindex;  // first CSR array
//    int* nlist;  // second CSR array
//    int* eweight;  // edge weights (-1, 0, 1)
//    int* origID;  // original node IDs
//};
//
//struct NodeData{
//    int nodeID;
//    int degree;
//    double status;
//    double authority;
//    double sumVacillation;
//    double influence;
//    int sumWeights;
//};
//
//struct EdgeData{
//    int destNodeID;
//    double vacillation;
//    double agreement;
//    double span;
//    double sumStatus;
//    double sumAuthority;
//};

static void freeGraph(Graph &g)
{
    g.nodes = 0;
    g.edges = 0;
    delete [] g.nindex;
    delete [] g.nlist;
    delete [] g.eweight;
    delete [] g.origID;
    g.nindex = NULL;
    g.nlist = NULL;
    g.eweight = NULL;
    g.origID = NULL;
}

static inline int representative(const int idx, int* const label)
{
    int curr = label[idx];
    if (curr != idx) {
        int next, prev = idx;
        while (curr > (next = label[curr])) {
            label[prev] = next;
            prev = curr;
            curr = next;
        }
    }
    return curr;
}

static Graph readGraph(const char* const name, std::vector<int> &nodeWeights)
{
    // read input from file
    FILE* fin = fopen(name, "rt");
    if (fin == NULL) {printf("ERROR: could not open input file %s\n", name); exit(-1);}
    size_t linesize = 256;
    char buf[linesize];
    char* ptr = buf;
    getline(&ptr, &linesize, fin);  // skip first line

    int selfedges = 0, wrongweights = 0, duplicates = 0, inconsistent = 0, line = 1, cnt = 0;
    int src, dst, wei;
    std::map<int, int> map;  // map node IDs to contiguous IDs
    std::set<std::pair<int, int>> set2;
    std::set<std::tuple<int, int, int>> set3;
    while (fscanf(fin, "%d,%d,%d", &src, &dst, &wei) == 3) {
        if (src == dst) {
            selfedges++;
        } else if ((wei < -1) || (wei > 1)) {  //MB
            wrongweights++;
        } else if (set2.find(std::make_pair(std::min(src, dst), std::max(src, dst))) != set2.end()) {
            if (set3.find(std::make_tuple(std::min(src, dst), std::max(src, dst), wei)) != set3.end()) {
                duplicates++;
            } else {
                inconsistent++;
            }
        } else {
            set2.insert(std::make_pair(std::min(src, dst), std::max(src, dst)));
            set3.insert(std::make_tuple(std::min(src, dst), std::max(src, dst), wei));
            if (map.find(src) == map.end()) {
                map[src] = cnt++;
            }
            if (map.find(dst) == map.end()) {
                map[dst] = cnt++;
            }
        }
        line++;
    }
    fclose(fin);

    std::set<std::tuple<int, int, int>>::iterator it;
    for(it = set3.begin(); it != set3.end(); it++){
        auto curNode = *it;
        int srcNode = std::get<0>(curNode);
        int dstNode = std::get<1>(curNode);
        int weight = std::get<2>(curNode);
        int max = std::max(srcNode, dstNode);

        if(max >= nodeWeights.size())
            nodeWeights.resize(max + 1);

        nodeWeights.at(srcNode) += weight;
        nodeWeights.at(dstNode) += weight;
    }

    // print stats
    printf("  read %d lines\n", line);
    if (selfedges > 0) printf("  skipped %d self-edges\n", selfedges);
    if (wrongweights > 0) printf("  skipped %d edges with out-of-range weights\n", wrongweights);
    if (duplicates > 0) printf("  skipped %d duplicate edges\n", duplicates);
    if (inconsistent > 0) printf("  skipped %d inconsistent edges\n", inconsistent);
    if (verify) {
        if ((int)map.size() != cnt) {printf("ERROR: wrong node count\n"); exit(-1);}
        printf("  number of unique nodes: %d\n", (int)map.size());
        printf("  number of unique edges: %d\n", (int)set3.size());
    }



    // compute CCs with union find
    int* const label = new int [cnt];
    for (int v = 0; v < cnt; v++) {
        label[v] = v;
    }
    for (auto ele: set3) {
        const int src = map[std::get<0>(ele)];
        const int dst = map[std::get<1>(ele)];
        const int vstat = representative(src, label);
        const int ostat = representative(dst, label);
        if (vstat != ostat) {
            if (vstat < ostat) {
                label[ostat] = vstat;
            } else {
                label[vstat] = ostat;
            }
        }
    }
    for (int v = 0; v < cnt; v++) {
        int next, vstat = label[v];
        while (vstat > (next = label[vstat])) {
            vstat = next;
        }
        label[v] = vstat;
    }

    // determine CC sizes
    int* const size = new int [cnt];
    for (int v = 0; v < cnt; v++) {
        size[v] = 0;
    }
    for (int v = 0; v < cnt; v++) {
        size[label[v]]++;
    }

    // find largest CC
    int hi = 0;
    for (int v = 1; v < cnt; v++) {
        if (size[hi] < size[v]) hi = v;
    }

    // keep if in largest CC and convert graph into set format
    Graph g;
    g.origID = new int [cnt];  // upper bound on size
    int nodes = 0, edges = 0;
    std::map<int, int> newmap;  // map node IDs to contiguous IDs
    std::set<std::pair<int, int>>* const node = new std::set<std::pair<int, int>> [cnt];  // upper bound on size
    for (auto ele: set3) {
        const int src = std::get<0>(ele);
        const int dst = std::get<1>(ele);
        const int wei = std::get<2>(ele);
        if (label[map[src]] == hi) {  // in largest CC
            if (newmap.find(src) == newmap.end()) {
                g.origID[nodes] = src;
                newmap[src] = nodes++;
            }
            if (newmap.find(dst) == newmap.end()) {
                g.origID[nodes] = dst;
                newmap[dst] = nodes++;
            }
            node[newmap[src]].insert(std::make_pair(newmap[dst], wei));
            node[newmap[dst]].insert(std::make_pair(newmap[src], wei));
            edges += 2;
        }
    }
    if (verify) {
        if (nodes > cnt) {printf("ERROR: too many nodes\n"); exit(-1);}
        if (edges > (int)set3.size() * 2) {printf("ERROR: too many edges\n"); exit(-1);}
    }

    // create graph in CSR format
    g.nodes = nodes;
    g.edges = edges;
    g.nindex = new int [g.nodes + 1];
    g.nlist = new int [g.edges];
    g.eweight = new int [g.edges];
    int acc = 0;
    for (int v = 0; v < g.nodes; v++) {
        g.nindex[v] = acc;
        for (auto ele: node[v]) {
            const int dst = ele.first;
            const int wei = ele.second;
            g.nlist[acc] = dst;
            g.eweight[acc] = wei;
            acc++;
        }
    }
    g.nindex[g.nodes] = acc;
    if (verify) {
        if (acc != edges) {printf("ERROR: wrong edge count in final graph\n"); exit(-1);}
    }

    delete [] label;
    delete [] size;
    delete [] node;

    return g;
}

// https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
static inline unsigned int hash(unsigned int val)
{
    val = ((val >> 16) ^ val) * 0x45d9f3b;
    val = ((val >> 16) ^ val) * 0x45d9f3b;
    return (val >> 16) ^ val;
}

static void init(const Graph& g, EdgeInfo* const einfo, int* const origEdge)
{
    // shift nlist, set minus if graph weight is -1, and init origEdge
#pragma omp parallel for default(none) shared(g)
    for (int j = 0; j < g.edges; j++) {
        g.nlist[j] <<= 1;
        einfo[j].end = (g.eweight[j] == -1) ? 1 : 0;
        origEdge[j] = j;
    }
}

static void generateSpanningTree(const Graph& g, const int root, const int seed, EdgeInfo* const einfo, int* const parent, int* const queue, int* const border, int* const label, int* const origEdge, bool* const edge_in_tree)
{
    const int seed2 = seed * seed + seed;

    // initialize
#pragma omp parallel for default(none) shared(g)
    for (int j = 0; j < g.edges; j++) g.nlist[j] &= ~1;  // edge is not in tree
#pragma omp parallel for default(none) shared(g)
    for (int i = 0; i < g.nodes; i++) parent[i] = -1;
    int tail = 1;
    parent[root] = INT_MAX & ~3;
    queue[0] = root;

    // BFS traversal
    int level = 0;
    border[0] = 0;
    border[1] = tail;
    while (border[level + 1] < g.nodes) {  // skip last iteration
        const int bit = (level & 1) | 2;
#pragma omp parallel for default(none) shared(g, level, tail)
        for (int i = border[level]; i < border[level + 1]; i++) {
            const int node = queue[i];
            const int me = (node << 2) | bit;
#pragma omp atomic write
            parent[node] = parent[node] & ~3;
            for (int j = g.nindex[node]; j < g.nindex[node + 1]; j++) {
                const int neighbor = g.nlist[j] >> 1;
                const int seed3 = neighbor ^ seed2;
                const int hash_me = hash(me ^ seed3);
                int val, hash_val;
                do {  // pick parent deterministically
#pragma omp atomic read
                    val = parent[neighbor];
                    hash_val = hash(val ^ seed3);
                } while (((val < 0) || (((val & 3) == bit) && ((hash_val < hash_me) || ((hash_val == hash_me) && (val < me))))) && (__sync_val_compare_and_swap(&parent[neighbor], val, me) != val));
                if (val < 0) {  // first visit
#pragma omp atomic capture
                    val = tail++;
                    queue[val] = neighbor;
                }
            }
        }
        level++;
        if (border[level] == tail) {printf("ERROR: input appears to have multiple connected components; terminating program\n"); exit(-1);}
        border[level + 1] = tail;
    }
    const int levels = level + 1;
    if (verify) {
        if (border[levels] != tail) {printf("ERROR: head mismatch\n"); exit(-1);}


        if (tail != g.nodes) {printf("ERROR: tail mismatch\n"); exit(-1);}
        for (int i = 0; i < g.nodes; i++) {
            if (parent[i] < 0) {printf("ERROR: found unvisited node %d\n", i); exit(-1);}
        }
    }

    // bottom up: push counts
#pragma omp parallel for default(none) shared(g)
    for (int i = 0; i < g.nodes; i++) label[i] = 1;
    for (int level = levels - 1; level > 0; level--) {  // skip level 0
#pragma omp parallel for default(none) shared(level)
        for (int i = border[level]; i < border[level + 1]; i++) {
            const int node = queue[i];
#pragma omp atomic
            label[parent[node] >> 2] += label[node];
        }
    }
    if (verify) {
        if (label[root] != g.nodes) {printf("ERROR: root count mismatch\n"); exit(-1);}
    }

    // top down: label tree + set nlist flag + set edge info + move tree nodes to front + make parent edge first in list
    label[root] = 0;
    for (int level = 0; level < levels; level++) {
#pragma omp parallel for default(none) shared(g, level)
        for (int i = border[level]; i < border[level + 1]; i++) {
            const int node = queue[i];
            const int par = parent[node] >> 2;
            const int nodelabel = label[node];
            const int beg = g.nindex[node];
            int paredge = -1;
            int lbl = (nodelabel >> 1) + 1;
            int pos = beg;
            for (int j = beg; j < g.nindex[node + 1]; j++) {
                const int neighbor = g.nlist[j] >> 1;
                if (neighbor == par) {
                    paredge = j;
                } else if ((parent[neighbor] >> 2) == node) {
                    const int count = label[neighbor];
                    label[neighbor] = lbl << 1;
                    lbl += count;
                    // set child edge info
                    einfo[j].beg = label[neighbor];
                    einfo[j].end = (einfo[j].end & 1) | ((lbl - 1) << 1);
                    g.nlist[j] |= 1;  // child edge is in tree
                    // swap
                    if (pos < j) {
                        std::swap(g.nlist[pos], g.nlist[j]);
                        std::swap(einfo[pos], einfo[j]);
                        std::swap(origEdge[pos], origEdge[j]);
                        if (paredge == pos) paredge = j;
                    }
                    pos++;
                }
            }
            if (paredge >= 0) {
                // set parent edge info
                einfo[paredge].beg = nodelabel | 1;
                einfo[paredge].end = (einfo[paredge].end & 1) | ((lbl - 1) << 1);
                g.nlist[paredge] |= 1;  // parent edge is in tree
                // move parent edge to front of list
                if (paredge != beg) {
                    if (paredge != pos) {
                        std::swap(g.nlist[pos], g.nlist[paredge]);
                        std::swap(einfo[pos], einfo[paredge]);
                        std::swap(origEdge[pos], origEdge[paredge]);
                        paredge = pos;
                    }
                    if (paredge != beg) {
                        std::swap(g.nlist[beg], g.nlist[paredge]);
                        std::swap(einfo[beg], einfo[paredge]);
                        std::swap(origEdge[beg], origEdge[paredge]);
                    }
                }
            }

            if (verify) {
                if (i == 0) {
                    if (lbl != g.nodes) {printf("ERROR: lbl mismatch\n"); exit(-1);}
                }
                int j = beg;
                while ((j < g.nindex[node + 1]) && (g.nlist[j] & 1)) j++;
                while ((j < g.nindex[node + 1]) && !(g.nlist[j] & 1)) j++;
                if (j != g.nindex[node + 1]) {printf("ERROR: not moved %d %d %d\n", beg, j, g.nindex[node + 1]); exit(-1);}
            }
        }
    }

    // update edge_in_tree
#pragma omp parallel for default(none) shared(g)
    for (int j = 0; j < g.edges; j++) {
        edge_in_tree[origEdge[j]] = g.nlist[j] & 1;
    }
}

static void initMinus(const Graph& g, const EdgeInfo* const einfo, bool* const minus)
{
    // set minus info to true
#pragma omp parallel for default(none) shared(g)
    for (int j = 0; j < g.edges; j++) {
        minus[j] = true;  //true so all duplicate edges will be excluded in determineCCs
    }

    // copy minus info of tree edges
#pragma omp parallel for default(none) shared(g) schedule(dynamic, 64)
    for (int i = 0; i < g.nodes; i++) {
        int j = g.nindex[i];
        while ((j < g.nindex[i + 1]) && (g.nlist[j] & 1)) {
            minus[j] = einfo[j].end & 1;
            j++;
        }
    }
}

static void processCycles(const Graph& g, const int* const label, EdgeInfo* const einfo, bool* const minus)
{
#pragma omp parallel for default(none) shared(g) schedule(dynamic, 64)
    for (int i = 0; i < g.nodes; i++) {
        const int target0 = label[i];
        const int target1 = target0 | 1;
        int j = g.nindex[i + 1] - 1;
        while ((j >= g.nindex[i]) && !(g.nlist[j] & 1)) {
            int curr = g.nlist[j] >> 1;
            if (curr > i) {  // only process edges in one direction
                int sum = 0;
                while (label[curr] != target0) {
                    int k = g.nindex[curr];
                    while ((einfo[k].beg & 1) == ((einfo[k].beg <= target1) && (target0 <= einfo[k].end))) k++;
                    if (verify) {
                        if ((k >= g.nindex[curr + 1]) || !(g.nlist[k] & 1)) {printf("ERROR: couldn't find path\n"); exit(-1);}
                    }
                    sum += einfo[k].end & 1;
                    curr = g.nlist[k] >> 1;
                }
                minus[j] = sum & 1;  // make cycle have even number of minuses
            }
            j--;
        }
    }
}

static void determineCCs(const Graph& g, int* const label, const bool* const minus, int* const count, const int* const origEdge, bool* const node_in_maj, bool* const node_has_neg_edge, bool* const edge_switched_sign)
{
    // init CCs
#pragma omp parallel for default(none) shared(g)
    for (int v = 0; v < g.nodes; v++) {
        label[v] = v;
        node_has_neg_edge[v] = false;
    }

    // compute CCs with union find
#pragma omp parallel for default(none) shared(g) schedule(dynamic, 64)
    for (int v = 0; v < g.nodes; v++) {
        const int beg = g.nindex[v];
        const int end = g.nindex[v + 1];
        int vstat = representative(v, label);
        for (int j = beg; j < end; j++) {
            const int nli = g.nlist[j] >> 1;
            if (v < nli) {
                edge_switched_sign[origEdge[j]] = ((minus[j] ? -1 : 1) != g.eweight[origEdge[j]]);  //MB: new
            }
            if (minus[j]) {
                if (v < nli) {
                    node_has_neg_edge[v] = true;
                    node_has_neg_edge[nli] = true;
                }
            } else {
                int ostat = representative(nli, label);
                bool repeat;
                do {
                    repeat = false;
                    if (vstat != ostat) {
                        int ret;
                        if (vstat < ostat) {
                            if ((ret = __sync_val_compare_and_swap(&label[ostat], ostat, vstat)) != ostat) {
                                ostat = ret;
                                repeat = true;
                            }
                        } else {
                            if ((ret = __sync_val_compare_and_swap(&label[vstat], vstat, ostat)) != vstat) {
                                vstat = ret;
                                repeat = true;
                            }
                        }
                    }
                } while (repeat);
            }
        }
    }

    // finalize CCs
#pragma omp parallel for default(none) shared(g)
    for (int v = 0; v < g.nodes; v++) {
        int next, vstat = label[v];
        const int old = vstat;
        while (vstat > (next = label[vstat])) {
            vstat = next;
        }
        if (old != vstat) label[v] = vstat;
    }

    // determine CC sizes
#pragma omp parallel for default(none) shared(g)
    for (int v = 0; v < g.nodes; v++) {
        count[v] = 0;
    }
#pragma omp parallel for default(none) shared(g)
    for (int v = 0; v < g.nodes; v++) {
#pragma omp atomic
        count[label[v]]++;
    }

    // find largest CC (source CC)
    int hi = 0;
#pragma omp parallel for default(none) shared(g, hi)
    for (int v = 1; v < g.nodes; v++) {
        if (count[hi] < count[v]) {
#pragma omp critical
            if (count[hi] < count[v]) hi = v;
        }
    }

    // init CC hop count (distance) from source CC, populate workset of edges that cross CCs
    std::set<std::pair<int, int>> ws;
    for (int v = 0; v < g.nodes; v++) {
        const int lblv = label[v];
        if (lblv == v) {
            count[lblv] = (lblv == hi) ? 0 : INT_MAX - 1;  // init count
        }
        for (int j = g.nindex[v]; j < g.nindex[v + 1]; j++) {
            const int nli = g.nlist[j] >> 1;
            const int lbln = label[nli];
            if (lblv < lbln) {  // only one direction
                ws.insert(std::make_pair(lblv, lbln));
            }
        }
    }

    // use Bellman Ford to compute distances
    bool changed;
    do {
        changed = false;
        for (auto p: ws) {
            const int lblv = p.first;
            const int lbln = p.second;
            const int distv = count[lblv];
            const int distn = count[lbln];
            if (distv + 1 < distn) {
                count[lbln] = distv + 1;
                changed = true;
            } else if (distn + 1 < distv) {
                count[lblv] = distn + 1;
                changed = true;
            }
        }
    } while (changed);

    // node is in majority if it is at even hop count from source CC
#pragma omp parallel for default(none) shared(g, hi)
    for (int v = 0; v < g.nodes; v++) {
        node_in_maj[v] = (count[label[v]] % 2) ^ 1;
    }
}

struct CPUTimer
{
    timeval beg, end;
    CPUTimer() {}
    ~CPUTimer() {}
    void start() {gettimeofday(&beg, NULL);}
    double elapsed() {gettimeofday(&end, NULL); return end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) / 1000000.0;}
};

int main(int argc, char* argv[])
{
    printf("Plus/Minus Code (%s)\n", __FILE__);
    printf("Copyright 2020 Texas State University\n");

    CPUTimer overall;
    overall.start();

    // process command line and read input
    if (argc != 3) {printf("USAGE: %s input_file_name iteration_count\n", argv[0]); exit(-1);}
    CPUTimer timer;
    timer.start();
    printf("verification: %s\n", verify ? "on" : "off");
    printf("input: %s\n", argv[1]);
    std::vector<int> nodeWeights(1);
    Graph g = readGraph(argv[1], nodeWeights);
    printf("nodes: %d\n", g.nodes);
    printf("edges: %d (%d)\n", g.edges, g.edges / 2);
    const int iterations = atoi(argv[2]);
    printf("trees: %d\n", iterations);
    printf("input time: %.6f s\n", timer.elapsed());

    // allocate all memory
    bool* const minus = new bool [g.edges];
    int* const parent = new int [g.nodes];
    int* const queue = new int [g.nodes];  // first used as queue, then as CC size
    int* const label = new int [g.nodes];  // first used as count, then as label, and finally as CC label
    int* const border = new int [g.nodes + 2];  // maybe make smaller
    int* const origEdge = new int [g.edges];  // original edge number (to undo shuffling)
    EdgeInfo* const einfo = new EdgeInfo [g.edges];
    int* const root = new int [g.nodes];

    // per tree stats
    bool* const node_in_maj = new bool [g.nodes * iterations]; //status
    bool* const node_has_neg_edge = new bool [g.nodes * iterations]; //authority
    bool* const edge_in_tree = new bool [g.edges * iterations]; //in spanning
    bool* const edge_switched_sign = new bool [g.edges * iterations]; //vacillation

    timer.start();
    for (int i = 0; i < g.nodes; i++) root[i] = i;
    if (pick_random_roots) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(root, root + g.nodes, gen);
    } else {
        std::partial_sort(root, root + std::min(iterations, g.nodes), root + g.nodes, [&](int a, int b) {
            return (g.nindex[a + 1] - g.nindex[a]) > (g.nindex[b + 1] - g.nindex[b]);
        });
    }
    init(g, einfo, origEdge);
    printf("init time:  %.6f s\n", timer.elapsed());

    for (int iter = 0; iter < iterations; iter++) {
        printf("tree %d\n", iter);

        // generate tree
        timer.start();
        generateSpanningTree(g, root[iter % g.nodes], iter + 17, einfo, parent, queue, border, label, origEdge, &edge_in_tree[iter * g.edges]);
        if (iter == 0) printf("tree time:  %.6f s\n", timer.elapsed());

        // initialize plus/minus
        timer.start();
        initMinus(g, einfo, minus);
        if (iter == 0) printf("pm time:    %.6f s\n", timer.elapsed());

        // find cycles
        timer.start();
        processCycles(g, label, einfo, minus);
        if (iter == 0) printf("cycle time: %.6f s\n", timer.elapsed());

        // determine connected components
        timer.start();
        determineCCs(g, label, minus, queue, origEdge, &node_in_maj[iter * g.nodes], &node_has_neg_edge[iter * g.nodes], &edge_switched_sign[iter * g.edges]);
        if (iter == 0) printf("CC time:    %.6f s\n", timer.elapsed());
    }

    // RESULT FILES SETUP
    std::string fullPath = argv[1];
    auto slashIndex = fullPath.rfind('/', fullPath.length());
    auto fileName = fullPath.substr(slashIndex + 1, fullPath.length() - slashIndex);
    auto uScoreIndex = fileName.rfind('_', fileName.length());
    auto outputFileName = fileName.substr(0, uScoreIndex);
    int numTrees = atoi(argv[2]);

    //STORING STATS FOR LATER
    int maxid = g.origID[0];
    for (int i = 1; i < g.nodes; i++) {
        maxid = std::max(maxid, g.origID[i]);
    }
    maxid++;
    std::vector<NodeData> nodeData(maxid);
    for (int i = 0; i < maxid; i++) {
        nodeData[i].nodeID = -1;
    }
    std::map<int, std::vector<EdgeData>> edgeData;


    //COMPUTE COLUMN FLAGS
#if DEBUG
    std::vector<int> onesLarger(iterations); //storing if edges w/ weight 1 are greater than edges w/ weight 0
#endif
    std::vector<bool> columnFlags(iterations);
    int numTies = 0;
    for(int j = 0; j < iterations; j++){
        int onesCount = 0;
        for(int k = 0; k < g.nodes; k++){
            onesCount += node_in_maj[j * g.nodes + k];
        }
        columnFlags[j] = (onesCount * 2 == g.nodes);
#if DEBUG
        if(onesCount * 2 > g.nodes){
            onesLarger[j] = 1;
        } else if(onesCount * 2 == g.nodes){
            onesLarger[j] = 0;
        } else{
            onesLarger[j] = -1;
        }
#endif
    }
#if DEBUG
    int numOnesLarger = 0;
    int numZerosLarger = 0;
    for(int i = 0; i < columnFlags.size(); i++){
        numTies += columnFlags[i];
        if(columnFlags[i] == 0 && onesLarger[i] == -1) {
            numZerosLarger++;
        } else if(columnFlags[i] == 0 && onesLarger[i] == 1)
            numOnesLarger++;
    }
    std::cout << "OnesLarger: " << numOnesLarger << " --- ZerosLarger: " << numZerosLarger << std::endl;
    std::cout << numTies << std::endl;
#endif
    //STATUS && AUTHORITY
    std::vector<int> nodesMajority(g.nodes);
    for(int i = 0; i < g.nodes; i++){
        int nodeStatus = 0;
        int nodeStatusSumAdd = 0;
        int nodeAuthority = 0;
        nodeData.at(g.origID[i]).nodeID = g.origID[i];
        for(int j = 0; j < iterations; j++){
            int curNode = node_in_maj[j * g.nodes + i];
            if(columnFlags[j])
                nodeStatus += 1;
            else
                nodeStatus += 2 * curNode;
            nodeStatusSumAdd += curNode;
            nodeAuthority += !node_has_neg_edge[j * g.nodes + i];
        }

        nodeData.at(g.origID[i]).degree = g.nindex[i + 1] - g.nindex[i];
        nodeData.at(g.origID[i]).status = nodeStatus / (double)(2 * numTrees);
        nodeData.at(g.origID[i]).authority = nodeAuthority / (double)numTrees;
    }

    for(int i = 0; i < g.nodes; i++)
        nodeData.at(g.origID[i]).sumWeights = nodeWeights.at(i);

#if TEST
    std::vector<std::pair<int, std::pair<int, int*>>> stateWeights; //weights of the created balanced states
    int* const iterWeights = new int[(iterations + 1) * g.edges];
    int weightPtrOffset = 0;
    std::vector<SpanEdge> spanningEdges; //spanning edges of the created balanced states
#endif

    // SPAN, VACILLATION, && AGREEMENT
    std::vector<int> treeWeightProducts(iterations);
    for (int v = 0; v < g.nodes; v++){
        for (int j = g.nindex[v]; j < g.nindex[v + 1]; j++){
            const int n = g.nlist[j] >> 1;
            if (v < n) {  // only print one copy of each edge (other copy does not hold correct data)
                int edgeSpan = 0;
                int edgeVac = 0;
                int edgeAgreement = 0;
#if TEST
                int curBaseWeight = g.eweight[origEdge[j]];
                stateWeights.emplace_back(std::make_pair(v, std::make_pair(n, nullptr)));
                stateWeights.back().second.second = iterWeights + weightPtrOffset;
                *(stateWeights.back().second.second) = curBaseWeight;
                weightPtrOffset++;
#endif
                for (int k = 0; k < iterations; k++) {
#if TEST
                    if(edge_in_tree[k * g.edges + origEdge[j]])
                        spanningEdges.emplace_back(k, v, n);
                    else
                        spanningEdges.emplace_back(-1, NULL, NULL);
                    *(iterWeights + weightPtrOffset++) = edge_switched_sign[k * g.edges + origEdge[j]] ? curBaseWeight * -1 : curBaseWeight;

#endif
                    edgeSpan += edge_in_tree[k * g.edges + origEdge[j]];
                    edgeVac += edge_switched_sign[k * g.edges + origEdge[j]];

                    bool srcNode = node_in_maj[k * g.nodes + v];
                    bool destNode = node_in_maj[k * g.nodes + n];
                    if(columnFlags[k])
                        edgeAgreement += 1;
                    else if(srcNode && destNode)
                        edgeAgreement += 2;

                    int treeWeight = 1;
                    for (int x = g.nindex[v]; x < g.nindex[v + 1]; x++) {
                        if(edge_in_tree[k * g.edges + origEdge[x]]){
                            int weight = g.eweight[origEdge[x]];
                            weight = (edge_switched_sign[k * g.edges + origEdge[x]]) ? weight * -1 : weight;
                            treeWeight *= weight;
                        }
                    }
                    treeWeightProducts.at(k) = treeWeight;
                }

                double normSpan = (double)edgeSpan / numTrees;
                double normVac = (double)edgeVac / numTrees;
                double normAgreement = (double)edgeAgreement / (2 * numTrees);
                auto foundEdgeSource = edgeData.find(g.origID[v]);
                double sumStatus = (nodeData.at(g.origID[n]).status + nodeData.at(g.origID[v]).status) / 2;
                EdgeData newConnection;
                newConnection.destNodeID = g.origID[n];
                newConnection.span = normSpan;
                newConnection.vacillation = normVac;
                newConnection.agreement = normAgreement;
                newConnection.sumStatus = sumStatus;

                if(foundEdgeSource == edgeData.end())
                    edgeData.insert({g.origID[v], std::vector<EdgeData>{newConnection}});
                else
                    foundEdgeSource->second.push_back(newConnection);
            }
        }
    }
    //put CSR printing here if needed from below

    //SUM VACILLATION AND INFLUENCE (SUM AGREEMENT)
    for(auto node: edgeData){
        for(auto connection: node.second){
            nodeData.at(node.first).sumVacillation += connection.vacillation;
            nodeData.at(connection.destNodeID).sumVacillation += connection.vacillation;
            nodeData.at(node.first).influence += connection.agreement;
            nodeData.at(connection.destNodeID).influence += connection.agreement;
        }
    }

    std::sort(stateWeights.begin(), stateWeights.end());
    for(auto edge: edgeData)
        std::sort(edge.second.begin(), edge.second.end(), [](const EdgeData& lhs, const EdgeData& rhs) { return lhs.destNodeID < rhs.destNodeID; });
    std::sort(spanningEdges.begin(), spanningEdges.end(), [](const SpanEdge& lhs, const SpanEdge& rhs){ return lhs.tree < rhs.tree; });

    std::string outputBasePath = R"(data-timing/Output_Data/Features/)" + outputFileName;
#if TEST
    std::ofstream nodeStat(R"(data-test\Output_Data\Features\)" + outputFileName + "\\" + outputFileName + "_nodeStat.csv");
    std::ofstream edgeStat(R"(data-test\Output_Data\Features\)" + outputFileName + "\\" + outputFileName +  "_edgeStat.csv");
    std::ofstream balancedStates(R"(data-test\Output_Data\Features\)" + outputFileName + "\\" + outputFileName +  "_balancedStates.csv");
    std::ofstream spanningTrees(R"(data-test\Output_Data\Features\)" + outputFileName + "\\" + outputFileName +  "_spanningTress.csv");
    std::ofstream oldNodePrints(R"(data-test\Output_Data\Features\)" + outputFileName + "\\" + outputFileName +  "_oldNodePrints.csv");
    std::ofstream oldEdgePrints(R"(data-test\Output_Data\Features\)" + outputFileName + "\\" + outputFileName +  "_oldEdgePrints.csv");
    std::ofstream treeWeightStat(R"(data-test\Output_Data\Features\)" + outputFileName + "\\" + outputFileName +  "_treeWeights.csv");
#else
    std::ofstream nodeStat(outputBasePath + "_nodeStat.csv");
    std::ofstream edgeStat(outputBasePath + "_edgeStat.csv");
    std::ofstream treeWeightStat(outputBasePath + "_treeWeights.csv");
#endif

    Printer::PrintNodeStats(nodeStat, nodeData, numTrees);
    Printer::PrintEdgeStats(edgeStat, edgeData, numTrees);
    Printer::PrintTreeWeights(treeWeightStat, treeWeightProducts);

#if TEST
    Printer::PrintBalanceStates(balancedStates, stateWeights, iterations);
    Printer::PrintSpanningTrees(spanningTrees, spanningEdges);
    Printer::PrintRawNodeArrays(oldNodePrints, g, node_in_maj, node_has_neg_edge, iterations);
    Printer::PrintRawEdgeArrays(oldEdgePrints, g, edge_in_tree, edge_switched_sign, origEdge, iterations);
#endif

    //PRINTOUTS FOR NINDEX, NLIST, EWEIGHTS
//    for(int i = 0; i < g.nodes + 1; i++)
//        std::cout << g.nindex[i] << " ";
//
//    std::cout << std::endl;
//    for(int i = 0 ; i < g.nodes; i++)
//    {
//        for(int j = g.nindex[i]; j < g.nindex[i + 1]; j++)
//            std::cout << (g.nlist[j] >> 1) << " ";
//    }
//    std::cout << std::endl;
//
//    for(int i = 0; i < g.nodes; i++)
//    {
//        for(int j = g.nindex[i]; j < g.nindex[i + 1]; j++)
//            std::cout << g.eweight[origEdge[j]] << " ";
//    }
//    std::cout << std::endl;
    //END PRINTOUTS


    // finalize
#if TEST
    balancedStates.close();
    oldEdgePrints.close();
    oldNodePrints.close();
#endif

    nodeStat.close();
    edgeStat.close();
    treeWeightStat.close();
    delete [] node_in_maj;
    delete [] node_has_neg_edge;
    delete [] edge_in_tree;
    delete [] edge_switched_sign;
    freeGraph(g);
    delete [] minus;
    delete [] einfo;
    delete [] parent;
    delete [] queue;
    delete [] label;
    delete [] border;
    delete [] origEdge;

    printf("overall runtime: %.6f s\n", overall.elapsed());
}

//data structure to hold largest CC data unwound from CSR format (mainly for testing)
//    std::vector<std::pair<int, std::vector<std::pair<int, int>>>> csrShift;
//    for (int v = 0; v < g.nodes; v++) {
//        std::cout << g.origID[v] << ": ";
//        auto newSrc = std::pair<int, std::vector<std::pair<int,int>>>();
//        newSrc.first = g.origID[v];
//        for (int j = g.nindex[v]; j < g.nindex[v + 1]; j++) {
//            const int n = g.nlist[j] >> 1;
//            auto newPairWeight = std::pair<int, int>(g.origID[n], g.eweight[origEdge[j]]);
//            newSrc.second.push_back(newPairWeight);
//        }
//        std::sort(newSrc.second.begin(), newSrc.second.end(), [](std::pair<int,int> fP, std::pair<int, int> sP){
//            return fP.first < sP.first;
//        });
//        csrShift.push_back(newSrc);
//    }
//    std::sort(csrShift.begin(), csrShift.end(), [](std::pair<int, std::vector<std::pair<int,int>>> fP, std::pair<int, std::vector<std::pair<int,int>>> sP){
//        return fP.first < sP.first;
//    });
//
//    //sum of nodeWeights
//    for(auto node: csrShift){
//        int sumNodeWeights = 0;
//        for(auto destNode: node.second){
//            sumNodeWeights += destNode.second;
//        }
//        nodeData.at(node.first).sumWeights = sumNodeWeights;
//    }
//
//    //csr printing area
//    for(auto node: csrShift){
//        std::cout << node.first << ": ";
//        for(auto destNode: node.second){
//            std::cout << "[" << destNode.first << ", " << destNode.second << "], ";
//        }
//        std::cout << std::endl;
//    }
