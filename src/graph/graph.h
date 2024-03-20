//
// Created by 高森森 on 2022/11/16.
//

#ifndef BASECINE_GRAPH_H
#define BASECINE_GRAPH_H

#include <queue>
#include <vector>
#include <string.h>
#include <iostream>
#include "../utils/types.h"
#include "../utils/utils.h"
#include "../utils/globals.h"
#include "MatchRecord.h"
#include "Neighbor.h"
#include "StarGraph.h"


class Edge{
private:
    uint v1;
    uint v2;
    uint v1Label;
    uint v2Label;
    uint eLabel;
    float eWeight;
    bool flag= true;
public:
    Edge(uint v1_, uint v2_, uint v1Label_, uint v2Label_, uint eLabel_,float eweight_):v1(v1_),v2(v2_),v1Label(v1Label_),v2Label(v2Label_),eLabel(eLabel_),eWeight(eweight_){
    }
    const uint GetV1() const{return this->v1;}
    const uint GetV2() const{return this->v2;}
    const uint GetV1Label() const{return this->v1Label;}
    const uint GetV2Label() const{return this->v2Label;}
    const uint GeteLabel() const{return this->eLabel;}
    const float GeteWeight()const{return this->eWeight;}
    bool operator == (const Edge& edge)const {
        if((edge.v1==this->v1&&edge.v2==this->v2&&edge.eLabel==this->eLabel)||
            (edge.v2==this->v1&&edge.v1==this->v2&&edge.eLabel==this->eLabel))
        {
            return true;
        }
        return false;
    }
    bool operator <(const Edge &edge)const{
        if(v1!=edge.v1){
            return v1<edge.v1;
        }
        else
            return v2<edge.v2;
    }
};





class Graph {
protected:
    uint edge_count_;//边个数
    uint vlabel_count_;//点标签个数
    uint elabel_count_;//边标签个数
    std::vector<std::vector<uint> > neighbors_;//邻接表
    std::vector<std::vector<uint> > elabels_;//边标签
    std::vector<std::vector<float>>weights_;//邻接表对应的边的权重
    std::vector<std::vector<uint>>timestamp_;//邻接表对应的边的时间戳
    std::vector<std::vector<vertexType>> matchVertexTypes;


public:
    std::queue<InsertUnit> updates_;//更新结构 插入或者删除边
    std::vector<uint> vlabels_;//节点标签
    std::vector<Edge>vEdge;//存储初始化载入所有的数据边
    std::vector<std::vector<Neighbor>> vNeighbors;//每个数据节点的邻居信息
    std::vector<int*>labelIndex;//记录每个节点邻居节点的点标签和边标签的个数
   // std::vector<std::vector<std::vector<int>>> LDRecord;
    //std::vector<std::vector<std::vector<ForwardNeighbor>>> forwardNeighbors;//前向邻居
    std::vector<std::vector<int>>isolatedRecord;//孤立节点的集合索引


public:
    Graph();

    virtual uint NumVertices() const { return vlabels_.size(); }
    virtual uint NumEdges() const { return edge_count_; }
    uint NumVLabels() const { return vlabel_count_; }
    uint NumELabels() const { return elabel_count_; }
    uint GetDiameter() const;

    void AddVertex(uint id, uint label);
    void RemoveVertex(uint id);
    void AddEdge(uint v1, uint v2, uint label,float weights,uint flag);//增加边的权重限制
    void RemoveEdge(uint flag,uint v1, uint v2);//flag=0表示是graphflow,flag=1 instopk

    uint GetVertexLabel(uint u) const;
    float GetEdgeWeight(uint v1,uint v2)const;
    uint GetEdgeTime(uint v1,uint v2)const;
    const std::vector<uint>& GetNeighbors(uint v) const;
    const std::vector<Neighbor>& GetGraphNeighbors(uint v) const;
    const std::vector<uint>& GetNeighborLabels(uint v) const;
    const std::vector<float>& GetNeighborWeights(uint v) const;
    uint GetDegree(uint v) const;
    std::tuple<uint, uint, uint> GetEdgeLabel(uint v1, uint v2) const;

    //0 query graph 1 data graph
    void LoadFromFile(const std::string &path,const uint flag);

    void createUpdateStream(const std::string data_path,const std::string dest_path);
    void createDataStream(const std::string data_path,const std::string dest_path);
    void createInitalDataGraph(const std::string data_path,const std::string dest_path);
    void LoadUpdateStream(const std::string &path);
    void PrintMetaData() const;
    void InitLabelIndex();
    void InitMatchOrderType( const std::vector<std::vector<uint> > &order_vs_,const std::vector<std::vector<std::vector<uint>>>&rightNeighbor);//初始化邻居标签分布
    void UpdateLabelIndex(uint v1,uint v2,uint label,uint flag);//flag=1为增加操作 flag=0为减少操作
    const vertexType GetVertexType(uint order_index,uint depth);
    std::vector<uint>GetIsolateVertexBeforeDepth(uint order_index,uint depth);
    void setBatchVertexType(uint order_index,const std::vector<uint> & vertexs, vertexType type);
    bool isNeighbor(uint v1,uint v2);
};


#endif //BASECINE_GRAPH_H
