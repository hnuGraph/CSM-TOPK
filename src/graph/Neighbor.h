//
// Created by ¸ß³þ³þ on 2023/3/23.
//

#ifndef BASELINE_NEIGHBOR_H
#define BASELINE_NEIGHBOR_H
#include "../utils/globals.h"
#include <utility>
class Neighbor{
protected:
    uint toVertexId;
    uint toVertexLabel;
    uint edgeLabel;
    float edgeWeight;
    uint matchQueryVertexId;
    uint fromQueryVertexId;
public:
    Neighbor();
    Neighbor(uint toVertexId_,uint toVertexLabel_,uint edgeLabel_,float edgeWeight_,uint matchQueryVertexId_,uint fromVertexId_):
    toVertexId(toVertexId_),
    toVertexLabel(toVertexLabel_),
    edgeLabel(edgeLabel_),
    edgeWeight(edgeWeight_),
    matchQueryVertexId(matchQueryVertexId_),
    fromQueryVertexId(fromVertexId_)
    {};
    Neighbor(uint toVertexId_,uint toVertexLabel_,float edgeLabel_):
    toVertexId(toVertexId_),
    toVertexLabel(toVertexLabel_),
    edgeLabel(edgeLabel_),
    edgeWeight(0){};
    Neighbor(uint toVertexId_,uint matchQueryVertexId_,uint fromVertexId_):
    toVertexId(toVertexId_),
    matchQueryVertexId(matchQueryVertexId_),
    fromQueryVertexId(fromVertexId_){};
    uint getVertexId() const;
    std::pair<uint,uint> GetelabelAndVertexLabel() const;
     bool operator>(const Neighbor &m) const;
    bool operator!=(const Neighbor&m)const;
     bool operator==(const Neighbor&m)const;
    float GetEdgeWeight()const;
    uint getVertexLabel()const;
    uint GetEdgelabel()const;
    const uint getMatchQueryVertexId()const ;
    const uint getfromVertexId()const ;

};

#endif //BASELINE_NEIGHBOR_H
