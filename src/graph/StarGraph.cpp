//
// Created by ¸ß³þ³þ on 2023/3/23.
//

#include "StarGraph.h"
void StarGraph::setStarMaxWeight(float w) {
    maxWeight=w;
}
const float StarGraph::getStarMaxWeight() {
    return maxWeight;
}
void StarGraph::InitalmaxWeight() {
    maxWeight=queryVertex.size()*mw;
}
const uint StarGraph::GetForwardNeighborNum() {
    return queryVertex.size();
}
void StarGraph::setMatchDataVertexId(uint id) {
    matchDataVertexId=id;
}
const uint StarGraph::getMatchDataVertexId() {
    return matchDataVertexId;
}
StarGraph::~StarGraph() {
    for(int i=0;i<queryVertex.size();i++){
        delete queryVertex[i];
        queryVertex[i]=NULL;
    }
    this->queryVertex.clear();
}
void StarGraph::computeMaxWeight() {
    maxWeight=0;
    for(auto f:queryVertex){
        maxWeight+=f->GetMaxWeight();
    }
}

std::vector<ForwardNeighbor *> &StarGraph::GetqueryVertex() {
return queryVertex;
}
