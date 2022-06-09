#include <set>
#include <climits>
#include <stack>
#include <queue>
#include <sys/time.h>
#include "ECLgraph.h"

// loops infinitely for inputs with multiple connected components

static const bool verification = false;

struct EdgeInfo {
  int beg;  // beginning of range
  int end;  // end of range
  bool inv;  // is range inverted or not
  signed char pm;  // plus or minus
  int plus;  // plus count
  int minus;  // minus count
};

static int* findBackEdges(const ECLgraph& g)
{
  // find opposite-direction edges
  int* const backEdge = new int [g.edges];
  if (verification) std::fill(backEdge, backEdge + g.edges, -1);

  for (int i = 0; i < g.nodes; i++) {
    for (int j = g.nindex[i]; j < g.nindex[i + 1]; j++) {
      const int neighbor = g.nlist[j];
      for (int k = g.nindex[neighbor]; k < g.nindex[neighbor + 1]; k++) {
        if (g.nlist[k] == i) {
          backEdge[j] = k;
          break;
        }
      }
    }
  }

  if (verification) {
    for (int j = 0; j < g.edges; j++) {
      if (backEdge[j] < 0) {printf("ERROR: missing back edge\n"); exit(-1);}
    }
  }

  return backEdge;
}

// source of hash function: https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
static unsigned int hash(unsigned int val)
{
  val = ((val >> 16) ^ val) * 0x45d9f3b;
  val = ((val >> 16) ^ val) * 0x45d9f3b;
  return (val >> 16) ^ val;
}

static bool* computeSpanningTree(const ECLgraph& g, const int* const backEdge, const int seed)
{
  const int root = hash(seed + 1) % g.nodes;

  // initialize
  bool* const treeEdge = new bool [g.edges];  // is edge in tree
  std::fill(treeEdge, treeEdge + g.edges, false);
  bool* const visited = new bool [g.nodes];
  std::fill(visited, visited + g.nodes, false);
  visited[root] = true;
  std::queue<int> q;
  q.push(root);

  // BFS traversal
  while (!q.empty()) {
    const int node = q.front();
    q.pop();
    for (int j = g.nindex[node + 1] - 1; j >= g.nindex[node]; j--) {  // reverse order
      const int neighbor = g.nlist[j];
      if (!visited[neighbor]) {
        visited[neighbor] = true;
        q.push(neighbor);
        treeEdge[j] = true;
        treeEdge[backEdge[j]] = true;
      }
    }
  }

  if (verification) {
    for (int i = 0; i < g.nodes; i++) {
      if (!visited[i]) printf("ERROR: found unvisited node %d\n", i);
    }
  }

  delete [] visited;
  return treeEdge;
}

static void moveTreeEdges(const ECLgraph& g, bool* const treeEdge, int* const backEdge, EdgeInfo* const einfo)
{
  for (int i = 0; i < g.nodes; i++) {
    int beg = g.nindex[i];
    int end = g.nindex[i + 1] - 1;
    while (beg < end) {
      while ((beg < end) && treeEdge[beg]) beg++;
      while ((beg < end) && !treeEdge[end]) end--;
      if (beg < end) {
        std::swap(treeEdge[beg], treeEdge[end]);
        std::swap(g.nlist[beg], g.nlist[end]);
        backEdge[backEdge[beg]] = end;
        backEdge[backEdge[end]] = beg;
        std::swap(backEdge[beg], backEdge[end]);
        std::swap(einfo[beg], einfo[end]);
        beg++;
        end--;
      }
    }

    if (verification) {
      int j = g.nindex[i];
      while ((j < g.nindex[i + 1]) && treeEdge[j]) j++;
      while ((j < g.nindex[i + 1]) && !treeEdge[j]) j++;
      if (j != g.nindex[i + 1]) {printf("ERROR: not moved %d %d %d\n", g.nindex[i], j, g.nindex[i + 1]); exit(-1);}
    }
  }
}

static int* labelTree(const ECLgraph& g, const bool* const treeEdge, const int* const backEdge, EdgeInfo* const einfo)
{
  int* const label = new int [g.nodes];  // new vertex ID
  std::fill(label, label + g.nodes, -1);

  // label tree in DFS preorder
  std::stack<int> edgestack;
  std::stack<bool> flagstack;
  int j = g.nindex[0];
  while ((j < g.nindex[1]) && treeEdge[j]) {
    edgestack.push(j);
    flagstack.push(false);
    j++;
  }
  label[0] = 0;
  int counter = 1;
  while (!edgestack.empty()) {
    const int edge = edgestack.top();
    const bool flag = flagstack.top();
    flagstack.pop();
    const int node = g.nlist[edge];
    if (flag) {
      edgestack.pop();
      einfo[edge].beg = label[node];
      einfo[edge].end = counter - 1;
      einfo[edge].inv = false;
      const int back = backEdge[edge];
      einfo[back].beg = label[node];//
      einfo[back].end = counter - 1;
      einfo[back].inv = true;
    } else {
      flagstack.push(true);
      label[node] = counter;
      counter++;
      int j = g.nindex[node];
      while ((j < g.nindex[node + 1]) && treeEdge[j]) {
        const int neighbor = g.nlist[j];
        if (label[neighbor] < 0) {
          edgestack.push(j);
          flagstack.push(false);
        }
        j++;
      }
    }
  }

  if (verification) {
    if (counter != g.nodes) printf("ERROR: counter = %d\n", counter);
    for (int i = 0; i < g.nodes; i++) {
      if (label[i] < 0) printf("ERROR: label[%d] = %d\n", i, label[i]);
    }
  }

  return label;
}

static void processCycles(const ECLgraph& g, const bool* const treeEdge, const int* const backEdge, const int* const label, EdgeInfo* const einfo)
{
  #pragma omp parallel for default(none) shared(g) schedule(dynamic)
  for (int i = 0; i < g.nodes; i++) {
    const int target = label[i];
    int j = g.nindex[i + 1] - 1;
    while ((j >= g.nindex[i]) && !treeEdge[j]) {
      int curr = g.nlist[j];
      if (curr > i) {  // only process edges in one direction
        int sum = 0;
        while (label[curr] != target) {
          const int oldcurr = curr;
          int k = g.nindex[curr];
          while (((einfo[k].inv || ((einfo[k].beg > target) || (target > einfo[k].end))) && (!einfo[k].inv || ((einfo[k].beg <= target) && (target <= einfo[k].end))))) k++;
          if (verification) {
            if ((k >= g.nindex[oldcurr + 1]) || !treeEdge[k]) {printf("ERROR: couldn't find path\n"); exit(-1);}
          }
          curr = g.nlist[k];
          sum += einfo[k].pm;
        }
        // some possibly wrong metric
        if (sum & 1) {
          einfo[j].plus++;
          einfo[backEdge[j]].plus++;
        } else {
          einfo[j].minus++;
          einfo[backEdge[j]].minus++;
        }
      }
      j--;
    }
  }

  if (verification) {
    for (int j = 0; j < g.edges; j++) {
      if (!treeEdge[j]) {
        if ((einfo[j].plus == 0) && (einfo[j].minus == 0)) {printf("ERROR: pm not set\n"); exit(-1);}
        if ((einfo[j].plus != einfo[backEdge[j]].plus) || (einfo[j].minus != einfo[backEdge[j]].minus)) {printf("ERROR: fwd and bwd do not match\n"); exit(-1);}
      }
    }
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
  CPUTimer overall;
  overall.start();

  if (argc != 3) {printf("USAGE: %s input_graph iteration_count\n", argv[0]); exit(-1);}
  printf("input: %s\n", argv[1]);
  ECLgraph g = readECLgraph(argv[1]);
  printf("nodes: %d\n", g.nodes);
  printf("edges: %d\n", g.edges);
  const int iterations = atoi(argv[2]);

  // compute back edges
  CPUTimer timer;
  timer.start();
  int* const backEdges = findBackEdges(g);
  printf("back time:  %.6f s\n", timer.elapsed());

  // use random plus minus
  timer.start();
  EdgeInfo* const einfo = new EdgeInfo [g.edges];
  for (int j = 0; j < g.edges; j++) {
    einfo[j].plus = 0;
    einfo[j].minus = 0;
    einfo[j].pm = hash(-1 - j) % 2;
    einfo[backEdges[j]].pm = einfo[j].pm;
  }
  printf("p/m time:   %.6f s\n", timer.elapsed());

  for (int iter = 0; iter < iterations; iter++) {
    printf("iteration %d\n", iter);

    // compute tree
    timer.start();
    bool* const treeEdges = computeSpanningTree(g, backEdges, iter);
    if (iter == 0) printf("tree time:  %.6f s\n", timer.elapsed());

    // move tree edges to front
    timer.start();
    moveTreeEdges(g, treeEdges, backEdges, einfo);
    if (iter == 0) printf("move time:  %.6f s\n", timer.elapsed());

    // label tree
    timer.start();
    const int* const label = labelTree(g, treeEdges, backEdges, einfo);
    if (iter == 0) printf("label time: %.6f s\n", timer.elapsed());

    // find cycles
    timer.start();
    processCycles(g, treeEdges, backEdges, label, einfo);
    if (iter == 0) printf("cycle time: %.6f s\n", timer.elapsed());

    delete [] treeEdges;
    delete [] label;
  }

  // print results
  for (int j = 0; j < 10/*g.edges*/; j++) {
    printf("%6d: %6d %6d\n", j, einfo[j].plus, einfo[j].minus);
  }

  freeECLgraph(g);
  delete [] backEdges;
  delete [] einfo;

  printf("overall runtime: %.6f s\n", overall.elapsed());
}

