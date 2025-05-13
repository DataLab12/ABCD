/* ABCD Version 2.0, 2024
 Author: Muhieddine Shebaro @ DataLab, Texas State University*/


#include <cstdio>
#include <climits>
#include <algorithm>
#include <set>
#include <map>
#include <sys/time.h>
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <queue>
#include <omp.h>
#include <random>

struct Graph {
  int nodes;
  int edges;
  int* nindex;  // first CSR array
  int* nlist;  // second CSR array
  int* eweight;  // edge weights (-1, 0, 1)
  int* origID;  // original node IDs
};

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

struct CPUTimer
{
  timeval beg, end;
  CPUTimer() {}
  ~CPUTimer() {}
  void start() {gettimeofday(&beg, NULL);}
  double elapsed() {gettimeofday(&end, NULL); return end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) / 1000000.0;}
};

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

static Graph readGraph_R(const char* const name,std::set<int> R)
{
  // read input from file
  FILE* fin = fopen(name, "rt");
  if (fin == NULL) {printf("ERROR: could not open input file %s\n", name); exit(-1);}
  size_t linesize = 256;
  char buf[linesize];
  char* ptr = buf;
  getline(&ptr, &linesize, fin);  // skip first line

  int selfedges = 0, wrongweights = 0, duplicates = 0, inconsistent = 0, line = 1, cnt = 0,removed_edges=0;
  int src, dst, wei;
  std::map<int, int> map;  // map node IDs to contiguous IDs
  std::set<std::pair<int, int>> set2;
  std::set<std::tuple<int, int, int>> set3;
  while (fscanf(fin, "%d,%d,%d", &src, &dst, &wei) == 3) {
    std::string total_string="";
          total_string.append(std::to_string(src));
          total_string.append(",");
          total_string.append(std::to_string(dst));
    if (src == dst) {
      selfedges++;
    } else if ((wei < -1) || (wei > 1)) {
      wrongweights++;
    }
    else if (R.find(int(src)) != R.end()|| R.find(int(dst)) != R.end()) {

           removed_edges++;
           continue;
                   }else if (set2.find(std::make_pair(std::min(src, dst), std::max(src, dst))) != set2.end()) {
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
  
  delete [] label;
  delete [] size;
  delete [] node;

  return g;
}

static Graph readGraph(const char* const name)
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
    } else if ((wei < -1) || (wei > 1)) {
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

  // print stats
  printf("  read %d lines\n", line);
  if (selfedges > 0) printf("  skipped %d self-edges\n", selfedges);
  if (wrongweights > 0) printf("  skipped %d edges with out-of-range weights\n", wrongweights);
  if (duplicates > 0) printf("  skipped %d duplicate edges\n", duplicates);
  if (inconsistent > 0) printf("  skipped %d inconsistent edges\n", inconsistent);
 

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
      if(wei==0){
        g.eweight[acc] = 1;

      }else{
        g.eweight[acc] = wei;

      }
      acc++;
    }
  }
    g.nindex[g.nodes] = acc;

  delete [] label;
  delete [] size;
  delete [] node;

  return g;
}


class Graphobj
{
public:
   int mSize;
    std::vector <int>* G;
    int n, m, t;
    //std::map<std::pair<int,int>,int> bridge;
    std::map<int,int> bridgeNode;
    std::vector<int> visited;
    std::vector<int> vk, l;
    std::vector<std::pair<int, int>> edges;

    Graphobj(int mSize,int n, int m) : mSize(mSize), n(n), m(m), visited(mSize), vk(mSize, -1), l(mSize, -1) {}




void depthSearch(int v, int p = -1) { 
   visited[v] = 1;
   vk[v] = l[v] = t++;
   for (auto x: G[v]) {
      if (x == p) {
         continue;
      }
      if (visited[x]) {
         l[v] = std::min(l[v], vk[x]);
      } else {
         depthSearch(x, v);
         l[v] = std::min(l[v], l[x]);
         if (l[x] > vk[v]) {
           //std::pair<int,int> p(x,v);
           //std::pair<int,int> p1(v,x);
            //   bridge[p] = 1;
           //    bridge[p1] = 1;
            bridgeNode[x]=1;
            bridgeNode[v]=1;
         }
      }
   }
}
void bridgeSearch() {
   t = 0;
   for (int i = 0; i < n; ++i) {
      if (!visited[i]) {
         depthSearch(i);
      }
   }
}
void solve(){
  G=new std::vector<int>[mSize];
   for (int i = 0; i < edges.size(); ++i) {
      int a, b; a = edges[i].first;
      b = edges[i].second;
      G[a].push_back(b);
      G[b].push_back(a);
   }
   bridgeSearch();
   delete [] G;
 
}
};




int main(int argc, char* argv[])
{
  printf("graphB++ balancing code for signed social network graphs (%s)\n", __FILE__);
  printf("Copyright 2024 Texas State University\n");

  CPUTimer overall;
  overall.start();
  Graph g;
  // process command line and read input
  if (argc != 3) {printf("USAGE: %s input_file_name iteration_count\n", argv[0]); exit(-1);}
  CPUTimer timer;
  timer.start();
  printf("input: %s\n", argv[1]);
  g=readGraph(argv[1]);
  printf("input time: %.6f s\n", timer.elapsed());
  int iterations=atoi(argv[2]);
  int iter=0;
  double lr=0.001;
    double* S=new double[g.nodes]; 
    
  
std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> distrib(-1, 1);
 
    for (int n = 0; n<g.nodes; n++){
      S[n]=distrib(gen);
       
    }

  int* const degree = new int [g.nodes];

    for (int v = 0; v < g.nodes; v++) {
       degree[v]=0;
     }

for (int v = 0; v < g.nodes; v++) {
    int index=g.nindex[v];
    int index_max=g.nindex[v+1];
    while(index<index_max){
      if(v<g.nlist[index]){
         degree[v]=degree[v]+1;
        degree[g.nlist[index]]=degree[g.nlist[index]]+1;
         
      }
    index=index+1;
    }
  }
      int frustration=0;

    while(iter<iterations){
        frustration=0;
        double* S_new=new double[g.nodes]; 
      for (int v = 0; v < g.nodes; v++) {
    const int beg = g.nindex[v];
    const int end = g.nindex[v + 1];
        double gradient=0.0;

    for (int j = beg; j < end; j++) {
      if(v<g.nlist[j]){
        int x,y;
        if(S[v]>=0){
          x=1;
        }else{
          x=-1;
        }
        if(S[g.nlist[j]]>=0){
          y=1;
        }else{
          y=-1;
        }
          frustration+=(1-x*y*g.eweight[j])/2;

      }

        gradient+=S[g.nlist[j]]*g.eweight[j];

    }
  gradient=-0.5*gradient;
    S_new[v]=S[v]-lr*gradient;

}
std::cout<<"Frustration: "<<frustration<<std::endl;
std::cout<<"Iteration: "<<iter+1<<std::endl;

for (int v = 0; v < g.nodes; v++) {
 S[v]=S_new[v];
}

delete [] S_new;


      iter++;
    }





Graphobj obj(g.nodes,g.nodes,g.edges/2);

 
         //std::cout<<obj.bridge[{3,2}]<<std::endl;
   
for (int v = 0; v < g.nodes; v++) {
    int index=g.nindex[v];
    int index_max=g.nindex[v+1];
    while(index<index_max){
      if(v<g.nlist[index]){
        obj.edges.push_back({v,g.nlist[index]});
      }
    index=index+1;
    }
  }
   obj.solve();


            
std::set<int> R;
int edges_handled=0;
for (int v = 0; v < g.nodes; v++) {
    const int beg = g.nindex[v];
    const int end = g.nindex[v + 1];
    for (int j = beg; j < end; j++) {
      if(v<g.nlist[j]){
        int x,y;
        if(S[v]>=0){
          x=1;
        }else{
          x=-1;
        }
        if(S[g.nlist[j]]>=0){
          y=1;
        }else{
          y=-1;
        }
          if((1-x*y*g.eweight[j])/2==1){

      if(R.find(g.origID[v])==R.end() && R.find(g.origID[g.nlist[j]])==R.end()){
  
          if(obj.bridgeNode[v]==1 && obj.bridgeNode[g.nlist[j]]==0){
            R.insert(g.origID[g.nlist[j]]);
            
          }
          if(obj.bridgeNode[v]==0 && obj.bridgeNode[g.nlist[j]]==1){
            R.insert(g.origID[v]);

          }
          if(obj.bridgeNode[v]==0 && obj.bridgeNode[g.nlist[j]]==0){
            if(degree[v]<degree[g.nlist[j]]){
            R.insert(g.origID[v]);
            }else{
              R.insert(g.origID[g.nlist[j]]);
            }

          }
          if(obj.bridgeNode[v]==1 && obj.bridgeNode[g.nlist[j]]==1){
            if(degree[v]<degree[g.nlist[j]]){
            R.insert(g.origID[v]);
            }else{
              R.insert(g.origID[g.nlist[j]]);
            }

          }
//std::cout<<obj.bridgeNode.size()<<std::endl;
/*for(auto const& x: obj.bridgeNode){
  std::cout<<x.first<<":"<<x.second<<std::endl;
}*/
/*obj.edges.clear();
obj.n=g.nodes;
obj.bridgeNode.clear();
obj.visited.clear();
obj.vk.clear();
obj.l.clear();
obj.vk.assign(obj.mSize, -1);
obj.l.assign(obj.mSize, -1);
obj.visited.assign(obj.mSize, 0);

for (int q = 0; q < g.nodes; q++) {
    int index=g.nindex[q];
    int index_max=g.nindex[q+1];
    while(index<index_max){
      if(q<g.nlist[index]){
if(R.find(g.origID[q]) == R.end()&& R.find(g.origID[g.nlist[index]]) == R.end()){
        obj.edges.push_back({q,g.nlist[index]});
}
      }
    index=index+1;
    }
  }

   obj.solve();
//std::cout<<obj.edges.size()<<std::endl;

    for (int q = 0; q < g.nodes; q++) {
       degree[q]=0;
     }

    for (int q = 0; q < g.nodes; q++) {
    int index=g.nindex[q];
    int index_max=g.nindex[q+1];
    while(index<index_max){
      if(q<g.nlist[index]){
        if(R.find(g.origID[q]) == R.end()&& R.find(g.origID[g.nlist[index]]) == R.end()){
         degree[q]=degree[q]+1;
        degree[g.nlist[index]]=degree[g.nlist[index]]+1;
        }
      }
    index=index+1;
    }
  }

*/

      }
edges_handled++;
std::cout<<"Edges handled: "<<edges_handled<<"/"<<frustration<<std::endl;



          }

      }


    }

}

Graph new_g=readGraph_R(argv[1],R);

 std::cout<<"Largest balanced sub-graph size: "<<new_g.nodes<<std::endl;

 freeGraph(new_g);
  delete [] S;
  delete [] degree;
  freeGraph(g);

  printf("overall runtime with I/O: %.6f s\n", overall.elapsed());
}
