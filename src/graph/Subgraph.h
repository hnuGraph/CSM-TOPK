//
// Created by ¸ß³þ³þ on 2023/6/18.
//

#ifndef BASELINE_SUBGRAPH_H
#define BASELINE_SUBGRAPH_H
#include "../utils/types.h"
#include "graph.h"
#include "unordered_set"

struct NeighborHash {
    std::size_t operator()(const Neighbor& neighbor) const {
        std::size_t seed = 0;
        seed ^= std::hash<uint>()(neighbor.getVertexId()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<uint>()(neighbor.getMatchQueryVertexId()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<uint>()(neighbor.getfromVertexId()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};
class Subgraph {
protected:
    uint edge_count_;

public:
    std::vector<std::vector<Neighbor>> vNeighbors;//data vertex neighbors
    std::vector<std::unordered_set<Neighbor,NeighborHash>>setVNeighbors;// Avoid duplicate neighboring vertexs
    std::vector<std::vector<uint>>matchCandidate;//Cadidate solutions for each node
    Subgraph(uint vertexSize,uint queryVertexSize);
    void addQueryVertexCandidate(uint q,uint v);
    bool AddEdge(uint u1,uint u2,uint v1,uint v1label,uint v2,uint v2label,uint label,float weight);
    bool RemoveEdge(uint v1,uint v1label, uint v2,uint v2label,uint u1,uint u2,uint elabel,uint weight);
    void deleteQueryVertexCandidate(uint q,uint v);
    const std::vector<Neighbor>& GetVNeighbors(uint v) const ;
};


#endif //BASELINE_SUBGRAPH_H
