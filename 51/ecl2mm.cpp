#include <stdlib.h>
#include <stdio.h>
#include "ECLgraph.h"

int main(int argc, char* argv[])
{
  printf("ECL to MatrixMarket Graph Converter (%s)\n", __FILE__);
  printf("Copyright 2016 Texas State University\n");

  if (argc != 3) {fprintf(stderr, "USAGE: %s input_file_name output_file_name\n\n", argv[0]);  exit(-1);}

  ECLgraph g = readECLgraph(argv[1]);

  printf("%s\t#name\n%d\t#nodes\n%d\t#edges\n", argv[1], g.nodes, g.edges);
  if (g.eweight != NULL) printf("yes\t#weights\n"); else printf("no\t#weights\n");

  FILE* fout = fopen(argv[2], "wt");
  if (fout == NULL) {fprintf(stderr, "ERROR: could not open output file %s\n\n", argv[2]);  exit(-1);}

  if (g.eweight != NULL) {
    fprintf(fout, "%%%%MatrixMarket matrix coordinate integer general\n");
  } else {
    fprintf(fout, "%%%%MatrixMarket matrix coordinate pattern general\n");
  }
  fprintf(fout, "%d\t%d\t%d\n", g.nodes, g.nodes, g.edges);

  for (int src = 0; src < g.nodes; src++) {
    for (int i = g.nindex[src]; i < g.nindex[src + 1]; i++) {
      int dst = g.nlist[i];
      if (g.eweight != NULL) {
        int wei = g.eweight[i];
        fprintf(fout, "%d\t%d\t%d\n", src + 1, dst + 1, wei);
      } else {
        fprintf(fout, "%d\t%d\n", src + 1, dst + 1);
      }
    }
  }
  fclose(fout);

  freeECLgraph(g);
  return 0;
}

