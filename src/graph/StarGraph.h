//
// Created by ¸ß³þ³þ on 2023/3/23.
//

#ifndef BASELINE_STARGRAPH_H
#define BASELINE_STARGRAPH_H
#include <vector>
#include "../utils/types.h"
#include "../utils/globals.h"
#include "ForwardNeighbor.h"


class StarGraph{
protected:
    std::vector<ForwardNeighbor*>queryVertex;
    float maxWeight;
    uint matchDataVertexId=UINT_MAX;
public:
    StarGraph(){};
    StarGraph(std::vector<ForwardNeighbor*>q){
        queryVertex=q;
    }
    ~StarGraph();
    void AddForwardNeighbor(ForwardNeighbor* f){
        queryVertex.emplace_back(f);
    }
    void InitalmaxWeight();
    void computeMaxWeight();
    const float getStarMaxWeight();
    void setStarMaxWeight(float w);
    void setMatchDataVertexId(uint id);
    const uint getMatchDataVertexId();
    const uint GetForwardNeighborNum();
    std::vector<ForwardNeighbor*>& GetqueryVertex();
};

#endif //BASELINE_STARGRAPH_H
