//
// Created by ¸ß³þ³þ on 2023/3/23.
//

#include "ForwardNeighbor.h"
const std::pair<uint,uint> ForwardNeighbor::GetelabelAndVertexLabel()  {
    return std::make_pair(edgeLabel,toVertexLabel);
}
const float ForwardNeighbor::GetMaxWeight() {
    return maxWeight;
}
 const uint ForwardNeighbor::GetVetexId() const{
    return toVertexId;
}
const uint ForwardNeighbor::GetElabel()const{
    return edgeLabel;
}
void ForwardNeighbor::SetMatchDataVertexId(uint id) {
    matchDataVertexId=id;
}
const uint ForwardNeighbor::GetMatchDataVertexId() const {
    return matchDataVertexId;
}
