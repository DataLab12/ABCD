#pragma once

#include <fstream>
#include <vector>
#include "DataStructures.h"
#include <iostream>
#include <map>

namespace Printer {

    void PrintNodeStats(std::ofstream& nodeStat, const std::vector<NodeData>& nodeData, int numTrees)
    {
        nodeStat << numTrees << " NodeID,Degree,Sum Weights,Status,Authority,Sum Vacillation,Sum Agreement\n";
        for(auto curNode: nodeData){
            nodeStat << curNode.nodeID << "," << curNode.degree << "," << curNode.sumWeights << "," << curNode.status << "," << curNode.authority
                     << "," << curNode.sumVacillation << ","<< curNode.influence << "\n";
        }
    }

    void PrintEdgeStats(std::ofstream& edgeStat, const std::map<int, std::vector<EdgeData>>& edgeData, int numTrees)
    {
        edgeStat << numTrees << " From NodeID,To NodeID,Span,Vacillation,Agreement,Sum Status\n";
        for(auto edge: edgeData){
            for(auto connection: edge.second){
                edgeStat << edge.first << "," << connection.destNodeID << "," << connection.span << "," <<
                         connection.vacillation << "," << connection.agreement << "," << connection.sumStatus << "\n";
            }
        }
    }

    void PrintTreeWeights(std::ofstream& treeWeightStat, const std::vector<int>& treeWeightProducts)
    {
        treeWeightStat << "Tree Number,Product of Tree Weights\n";
        for(int i = 0; i < treeWeightProducts.size(); i++)
            treeWeightStat << i << "," << treeWeightProducts.at(i) << "\n";
    }

    void PrintBalanceStates(std::ofstream& balancedStates, const std::vector<std::pair<int, std::pair<int, int*>>>& stateWeights, int iterations)
    {
        balancedStates << "Edge Src,Edge Dest,Beginning State";
        for(int i = 0; i < iterations; i++)
            balancedStates << ",S" << i;
        balancedStates << "\n";

        for(auto node: stateWeights)
        {
            balancedStates << node.first << "," << node.second.first << "," << *(node.second.second);
            for(int i = 0; i < iterations; i++)
                balancedStates << "," << *(node.second.second + i + 1);
            balancedStates << "\n";
        }
    }

    void PrintSpanningTrees(std::ofstream &spanningTrees, std::vector<SpanEdge> spanningEdges)
    {
        spanningTrees << "Iteration,Path\n";
        int curTree;
        int i = 0;
        while(spanningEdges.at(i).tree == -1){ i++; }
        for(;i < spanningEdges.size();)
        {
            curTree = spanningEdges.at(i).tree;
            spanningTrees << curTree;
            while(i < spanningEdges.size() && spanningEdges.at(i).tree == curTree)
            {
                spanningTrees << "," << spanningEdges.at(i).src << "->" << spanningEdges.at(i).dest;
                i++;
            }
            spanningTrees << "\n";
        }

    }

    void PrintRawNodeArrays(std::ofstream& oldNodePrints, const Graph& g, const bool* node_in_maj, const bool* node_has_neg_edge, int iterations)
    {
        oldNodePrints << "original node ID, node_in_maj [for each tree]\n";
        for (int i = 0; i < g.nodes; i++)
        {
            oldNodePrints << g.origID[i];
            for (int j = 0; j < iterations; j++)
                oldNodePrints << node_in_maj[j * g.nodes + i];
            oldNodePrints << "\n";
        }

        oldNodePrints << "original node ID, node_has_neg_edge [for each tree]\n";
        for (int i = 0; i < g.nodes; i++)
        {
            oldNodePrints << g.origID[i];
            for (int j = 0; j < iterations; j++)
                oldNodePrints << node_has_neg_edge[j * g.nodes + i];
            oldNodePrints << "\n";
        }
    }

    void PrintRawEdgeArrays(std::ofstream& oldEdgePrints, const Graph& g, const bool* edge_in_tree, const bool* edge_switched_sign, const int* origEdge, int iterations)
    {
        oldEdgePrints << "source node ID, destination node ID, edge_in_tree [for each tree]\n";
        for (int v = 0; v < g.nodes; v++)
        {
            for (int j = g.nindex[v]; j < g.nindex[v + 1]; j++) {
                const int n = g.nlist[j] >> 1;
                if (v < n) {  // only print one copy of each edge (other copy does not hold correct data)
                    oldEdgePrints << g.origID[v] << "," << g.origID[n];
                    for (int k = 0; k < iterations; k++) {
                        oldEdgePrints << "," << edge_in_tree[k * g.edges + origEdge[j]];
                    }
                    oldEdgePrints << "\n";
                }
            }
        }

        oldEdgePrints << "source node ID, destination node ID, edge_switched_sign [for each tree]\n";
        for (int v = 0; v < g.nodes; v++)
        {
            for (int j = g.nindex[v]; j < g.nindex[v + 1]; j++) {
                const int n = g.nlist[j] >> 1;
                if (v < n) {  // only print one copy of each edge (other copy does not hold correct data)
                    oldEdgePrints << g.origID[v] << "," << g.origID[n];
                    for (int k = 0; k < iterations; k++) {
                        oldEdgePrints << "," << edge_switched_sign[k * g.edges + origEdge[j]];
                    }
                    oldEdgePrints << "\n";
                }
            }
        }
    }

}