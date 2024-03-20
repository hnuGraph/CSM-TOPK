//
// Created by 高楚楚 on 2023/6/18.
//

#include "Subgraph.h"



Subgraph::Subgraph(uint vertexSize,uint queryVertexSize):
vNeighbors(vertexSize),
matchCandidate(queryVertexSize),
setVNeighbors(vertexSize)
{}

void Subgraph::addQueryVertexCandidate(uint q, uint v) {
    //按序插入
    auto it = std::lower_bound(matchCandidate[q].begin(), matchCandidate[q].end(), v);
    // 在确定的位置插入新的元素
    matchCandidate[q].insert(it, v);
//    matchCandidate[q].emplace_back(v);
}
bool Subgraph::AddEdge(uint u1,uint u2,uint v1, uint v1label,uint v2,uint v2label, uint label, float weight) {
//找到大于等于v2的迭代器的位置
    Neighbor neighbor(v2,v2label,label,weight,u2,u1);
    if(setVNeighbors[v1].count(neighbor)){
        return false;
    }
    else{
        auto lower=std::lower_bound(vNeighbors[v1].begin(),vNeighbors[v1].end(),neighbor,std::greater<Neighbor>());
        vNeighbors[v1].insert(lower,neighbor);
        //vNeighbors[v1].emplace_back(neighbor);
        setVNeighbors[v1].insert(neighbor);
    }
    Neighbor neighbor2(v1, v1label,label,weight,u1,u2);
    if(setVNeighbors[v2].count(neighbor2)){
        return false;
    }
    else{
        auto lower=std::lower_bound(vNeighbors[v2].begin(),vNeighbors[v2].end(),neighbor2,std::greater<Neighbor>());
        vNeighbors[v2].insert(lower,neighbor2);
        //Neighbors[v2].emplace_back(neighbor2);
        setVNeighbors[v2].insert(neighbor2);
    }
    edge_count_++;
    return true;
}

bool Subgraph::RemoveEdge(uint v1,uint v1label, uint v2,uint v2label,uint u1,uint u2,uint elabel,uint weight)
{
    Neighbor neighbor(v2,v2label,elabel,weight,u2,u1);
    if(!setVNeighbors[v1].count(neighbor)){
        return false;
    }
    auto lower=std::lower_bound(vNeighbors[v1].begin(),vNeighbors[v1].end(),neighbor,std::greater<Neighbor>());
    if(*lower!=neighbor){
        std::cout<<"deletion error"<<std::endl;
    }
    vNeighbors[v1].erase(lower);
    /*for(auto it=vNeighbors[v1].begin();it!=vNeighbors[v1].end();it++){
        if((*it)==neighbor){
            vNeighbors[v1].erase(it);
            break;
        }
    }*/
    setVNeighbors[v1].erase(neighbor);

    Neighbor neighbor2(v1,v1label,elabel,weight,u1,u2);
    if(!setVNeighbors[v2].count(neighbor2)){
        return false;
    }
    auto lower2=std::lower_bound(vNeighbors[v2].begin(),vNeighbors[v2].end(),neighbor2,std::greater<Neighbor>());
    if(lower2 == vNeighbors[v2].end()||*lower2!=neighbor2){
        std::cout<<"deletion error"<<std::endl;
        exit(-1);
    }
    vNeighbors[v2].erase(lower2);
  /*  for(auto it=vNeighbors[v2].begin();it!=vNeighbors[v2].end();it++){
        if((*it)==neighbor2){
            vNeighbors[v2].erase(it);
            break;
        }
    }*/
    setVNeighbors[v2].erase(neighbor2);
    edge_count_--;
    return true;
}
void Subgraph::deleteQueryVertexCandidate(uint q, uint v) {
    std::vector<uint>&candidate=matchCandidate[q];
    auto it=std::lower_bound(candidate.begin(),candidate.end(),v);

    if (it != candidate.end() && *it ==v) {
        // 删除目标元素
        candidate.erase(it);
    }


}
const std::vector<Neighbor>& Subgraph::GetVNeighbors(uint v) const {
    return vNeighbors[v];
}



