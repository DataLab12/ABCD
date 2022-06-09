#include <vector>
#include <map>
#include <set>
#include "ECLgraph.h"

int main(int argc, char* argv[])
{
  printf("CSV to ECL Graph Converter (%s)\n", __FILE__);
  printf("Copyright 2020 Texas State University\n");

  if (argc != 3) {fprintf(stderr, "USAGE: %s input_file_name output_file_name\n\n", argv[0]);  exit(-1);}

  // read file into vector
  FILE* fin = fopen(argv[1], "rt");  if (fin == NULL) {fprintf(stderr, "ERROR: could not open input file %s\n\n", argv[1]);  exit(-1);}
  char line[256];
  char* ptr = line;
  size_t linesize = 256;
  getline(&ptr, &linesize, fin);
  int src, dst, wei;
  std::vector<std::tuple<int, int, int>> vec;
  while (fscanf(fin, "%d,%d,%d", &src, &dst, &wei) == 3) {
    vec.push_back(std::make_tuple(src, dst, wei));
  }
  fclose(fin);
  printf("number of edges: %d\n", (int)vec.size() * 2);

  // map node IDs to contiguous IDs
  int cnt = 0;
  std::map<int, int> map;
  for (auto ele: vec) {
    const int src = std::get<0>(ele);
    const int dst = std::get<1>(ele);
    if (map.find(src) == map.end()) map[src] = cnt++;
    if (map.find(dst) == map.end()) map[dst] = cnt++;
  }
  printf("number of nodes: %d %d\n", (int)map.size(), cnt);

  // create graph in set format
  std::set<std::pair<int, int>>* const node = new std::set<std::pair<int, int>> [cnt];
  for (auto ele: vec) {
    const int src = std::get<0>(ele);
    const int dst = std::get<1>(ele);
    const int wei = std::get<2>(ele);
    node[map[src]].insert(std::make_pair(map[dst], wei));
    node[map[dst]].insert(std::make_pair(map[src], wei));
  }

  // create graph in ECL format
  ECLgraph g;
  g.nodes = cnt;
  g.edges = 2 * vec.size();
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
  printf("edges filled: %d\n", acc);
  g.edges = acc;  // in case of duplicates

  writeECLgraph(g, argv[2]);
  freeECLgraph(g);
  delete [] node;

  return 0;
}
