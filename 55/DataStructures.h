//
// Created by brhod on 2/14/2021.
//

#ifndef STATUS_DATASTRUCTURES_H
#define STATUS_DATASTRUCTURES_H

struct Graph {
    int nodes;
    int edges;
    int* nindex;  // first CSR array
    int* nlist;  // second CSR array
    int* eweight;  // edge weights (-1, 0, 1)
    int* origID;  // original node IDs
};

struct NodeData{
    int nodeID;
    int degree;
    double status;
    double authority;
    double sumVacillation;
    double influence;
    int sumWeights;
};

struct EdgeData{
    int destNodeID;
    double vacillation;
    double agreement;
    double span;
    double sumStatus;
    double sumAuthority;
};

struct SpanEdge{
    int tree;
    int src;
    int dest;
    SpanEdge(int t, int s, int d) : tree(t), src(s), dest(d){}
};

#endif //STATUS_DATASTRUCTURES_H
