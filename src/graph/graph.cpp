#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <vector>
#include "graph.h"
#include "../utils/types.h"
#include "../utils/utils.h"
#include "../utils/globals.h"

Graph::Graph()
: edge_count_(0)
, vlabel_count_(0)
, elabel_count_(0)
, neighbors_{}
, elabels_{}
, updates_{}
, vlabels_{}
{}

void Graph::AddVertex(uint id, uint label)
{
    //节点数据必须是有序的
    if (id >= vlabels_.size())
    {
        vlabels_.resize(id + 1, NOT_EXIST);
        vlabels_[id] = label;
        neighbors_.resize(id + 1);
        elabels_.resize(id + 1);
        weights_.resize(id+1);
        timestamp_.resize(id+1);
        vNeighbors.resize(id+1);
    }
    else if (vlabels_[id] == NOT_EXIST)
    {
        vlabels_[id] = label;
    }
    //保证节点的标签按顺序插入
    vlabel_count_ = std::max(vlabel_count_, label + 1);
    // print graph
    /*std::cout << "labels: ";
    for (uint i = 0; i < vlabels_.size(); i++)
    {
        std::cout << i << ":" << vlabels_[i] << " (";
        for (uint j = 0; j < neighbors_[i].size(); j++)
        {
            std::cout << neighbors_[i][j] << ":" << elabels_[i][j] << " ";
        }
        std::cout << ")" << std::endl;
    }*/
}

void Graph::RemoveVertex(uint id)
{
    vlabels_[id] = NOT_EXIST;
    neighbors_[id].clear();
    elabels_[id].clear();
    weights_[id].clear();
    timestamp_[id].clear();
    vNeighbors[id].clear();
}

void Graph::AddEdge(uint v1, uint v2, uint label,float weights,uint flag)
{
    //找到大于等于v2的迭代器的位置
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower != neighbors_[v1].end() && *lower == v2) return;
    
    size_t dis = std::distance(neighbors_[v1].begin(), lower);
    neighbors_[v1].insert(lower, v2);
    elabels_[v1].insert(elabels_[v1].begin() + dis, label);
    Neighbor neighbor(v2,this->GetVertexLabel(v2),label,weights,0,0);
    vNeighbors[v1].emplace_back(neighbor);
    if(flag)
    {
        weights_[v1].insert(weights_[v1].begin()+dis,weights);
    }

    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    if (lower != neighbors_[v2].end() && *lower == v1) return;
    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].insert(lower, v1);
    elabels_[v2].insert(elabels_[v2].begin() + dis, label);
    Neighbor neighbor2(v1, this->GetVertexLabel(v1),label,weights,0,0);
    vNeighbors[v2].emplace_back(neighbor2);
    if(flag)
    {
        weights_[v2].insert(weights_[v2].begin()+dis,weights);
    }


    if(v1>v2)
    {
        std::swap(v1,v2);
    }
    Edge edge(v1,v2,this->GetVertexLabel(v1),this->GetVertexLabel(v2),label,weights);
    this->vEdge.emplace_back(edge);
    edge_count_++;
    elabel_count_ = std::max(elabel_count_, label + 1);
    // print graph
    /*std::cout << "labels: ";
    for (uint i = 0; i < vlabels_.size(); i++)
    {
        std::cout << i << ":" << vlabels_[i] << " (";
        for (uint j = 0; j < neighbors_[i].size(); j++)
        {
            std::cout << neighbors_[i][j] << ":" << elabels_[i][j] << " ";
        }
        std::cout << ")" << std::endl;
    }*/
}

void Graph::RemoveEdge(uint flag,uint v1, uint v2)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower == neighbors_[v1].end() || *lower != v2)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v1].erase(lower);
    elabels_[v1].erase(elabels_[v1].begin() + std::distance(neighbors_[v1].begin(), lower));
    weights_[v1].erase(weights_[v1].begin()+std::distance(neighbors_[v1].begin(),lower));
    vNeighbors[v1].erase(vNeighbors[v1].begin()+std::distance(neighbors_[v1].begin(),lower));
   /* for(auto it=vNeighbors[v1].begin();it!=vNeighbors[v1].end();it++){
        if(it->getVertexId()==v2){
           vNeighbors[v1].erase(it);
            break;
        }
    }*/
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    if (lower == neighbors_[v2].end() || *lower != v1)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v2].erase(lower);
    elabels_[v2].erase(elabels_[v2].begin() + std::distance(neighbors_[v2].begin(), lower));
    weights_[v2].erase(weights_[v2].begin()+std::distance(neighbors_[v2].begin(),lower));
    vNeighbors[v2].erase(vNeighbors[v2].begin()+std::distance(neighbors_[v2].begin(),lower));
    if(flag){
        if(v1>v2){
            std::swap(v1,v2);
        }
        for(int i=0;i<vEdge.size();i++){
            const Edge &edge=vEdge[i];
            if(edge.GetV1()==v1&&edge.GetV2()==v2){
                vEdge.erase(vEdge.begin()+i);
                break;
            }
        }
    }
    edge_count_--;
}

uint Graph::GetVertexLabel(uint u) const
{
    return vlabels_[u];
}

const std::vector<uint>& Graph::GetNeighbors(uint v) const
{
    return neighbors_[v];
}
const std::vector<Neighbor>& Graph::GetGraphNeighbors(uint v) const {
    return vNeighbors[v];
}
const std::vector<uint>& Graph::GetNeighborLabels(uint v) const
{
    return elabels_[v];
}
const std::vector<float>& Graph::GetNeighborWeights(uint v) const {
    return weights_[v];
}
float Graph::GetEdgeWeight(uint v1,uint v2) const {
    float e_weight;

    const std::vector<uint> *nbrs;
    const std::vector<float> *weights;
    uint other;
    //从邻居少的节点来找边标签
    if (GetDegree(v1) < GetDegree(v2))
    {
        nbrs = &GetNeighbors(v1);
        weights = &weights_[v1];
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        weights= &weights_[v2];
        other = v1;
    }

    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)
    {
        mid = (start + end) / 2;
        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            e_weight = weights->at(mid);
            return e_weight;
        }
    }
    return 0;
}

uint Graph::GetEdgeTime(uint v1, uint v2) const {
    uint e_time;

    const std::vector<uint> *nbrs;
    const std::vector<uint> *times;
    uint other;
    //从邻居少的节点来找边标签
    if (GetDegree(v1) < GetDegree(v2))
    {
        nbrs = &GetNeighbors(v1);
        times = &timestamp_[v1];
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        times= &timestamp_[v2];
        other = v1;
    }

    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)
    {
        mid = (start + end) / 2;
        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            e_time = times->at(mid);
            return e_time;
        }
    }
    return 0;
}

std::tuple<uint, uint, uint> Graph::GetEdgeLabel(uint v1, uint v2) const
{
    uint v1_label, v2_label, e_label;
    v1_label = GetVertexLabel(v1);
    v2_label = GetVertexLabel(v2);

    const std::vector<uint> *nbrs;
    const std::vector<uint> *elabel;
    uint other;
    //从邻居少的节点来找边标签
    if (GetDegree(v1) < GetDegree(v2))
    {
        nbrs = &GetNeighbors(v1);
        elabel = &elabels_[v1];
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        elabel = &elabels_[v2];
        other = v1;
    }
    
    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)
    {
        mid = (start + end) / 2;
        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            e_label = elabel->at(mid);
            return std::tuple<uint,uint,uint>{v1_label, v2_label, e_label};
        }
    }
    return std::tuple<uint,uint,uint>{v1_label, v2_label, -1};
}

uint Graph::GetDegree(uint v) const
{
    return neighbors_[v].size();
}

uint Graph::GetDiameter() const
{
    uint diameter = 0;
    for (uint i = 0u; i < NumVertices(); i++)
    if (GetVertexLabel(i) != NOT_EXIST)
    {
        std::queue<uint> bfs_queue;
        std::vector<bool> visited(NumVertices(), false);
        uint level = UINT_MAX;
        bfs_queue.push(i);
        visited[i] = true;
        while (!bfs_queue.empty())
        {
            level++;
            uint size = bfs_queue.size();
            for (uint j = 0u; j < size; j++)
            {
                uint front = bfs_queue.front();
                bfs_queue.pop();

                const auto& nbrs = GetNeighbors(front);
                for (const uint nbr: nbrs)
                {
                    if (!visited[nbr])
                    {
                        bfs_queue.push(nbr);
                        visited[nbr] = true;
                    }
                }
            }
        }
        if (level > diameter) diameter = level;
    }
    return diameter;
}
//待改，将查询图和数据图得分开
void Graph::LoadFromFile(const std::string &path,const uint flag)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    char type;
    while (ifs >> type)
    {
        if (type == 't')
        {
            char temp1;
            uint temp2;
            ifs >> temp1 >> temp2;
        }
        else if (type == 'v')
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            AddVertex(vertex_id, label);
        }
        else
        {
            uint from_id, to_id, label,timestamp;
            float weights=0;
            timestamp=0;
            label=0;
            if(flag==0){
                ifs>>from_id>>to_id>>label;
                uint f=flag;
                AddEdge(from_id,to_id,label,weights,f);
            }else{
                ifs >> from_id >> to_id >> label>>weights;
                AddEdge(from_id, to_id, label,weights,flag);
            }

        }
    }
    ifs.close();
}


void Graph::createDataStream(const std::string data_path, const std::string dest_path) {
    std::ifstream ifs(data_path);
    std::ofstream ofs(dest_path,std::ios::out);
    if (!ifs.is_open())
    {
        std::cout << "Failed to open: " <<data_path << std::endl;
        exit(-1);
    }
    std::string line;
    float weight=0;
    while(std::getline(ifs,line)){
        weight=rand()%100;
        if(line[0]=='v'){
            ofs<<line<<std::endl;

        }else{
//            line.erase(line.length()-1);
            ofs<<line<<" "<<weight<<" "<<ts<<std::endl;
            ts++;
        }
    }
    ifs.close();
    ofs.close();
}
void Graph::createInitalDataGraph(const std::string data_path, const std::string dest_path) {
    std::ifstream ifs(data_path);
    std::ofstream ofs(dest_path,std::ios::out);
    if (!ifs.is_open())
    {
        std::cout << "Failed to open: " <<data_path << std::endl;
        exit(-1);
    }
    std::string line;
    uint cnt=0;
    while(cnt<windowSize){
        std::getline(ifs,line);
        if(line[0]=='v'){
            ofs<<line<<std::endl;
        }else{
            ofs<<line<<std::endl;
            cnt++;
        }
    }
    ifs.close();
    ofs.close();
}
void Graph::createUpdateStream(const std::string data_path, const std::string dest_path) {
    std::ifstream ifs(data_path);
    std::ofstream ofs(dest_path,std::ios::out);
    if (!ifs.is_open())
    {
        std::cout << "Failed to open: " <<data_path << std::endl;
        exit(-1);
    }
    std::string line;
    int cnt=0;
    std::vector<std::string>lines;
    while(std::getline(ifs,line)){
        lines.emplace_back(line);
    }
    auto it1=lines.begin();
    auto it2=lines.end();
    for(auto it=lines.begin();it!=lines.end();it++){
        if(it->c_str()[0]=='e'){
            it1=it;
            it2=it+windowSize;
            break;
        }
    }
/*    std::cout<<it1->c_str()<<std::endl;
    std::cout<<it2->c_str()<<std::endl;*/

    while(it2!=lines.end()){
        ofs<<it2->c_str()<<std::endl;
        ofs<<"-"<<it1->c_str()<<std::endl;
        it1++;
        it2++;
    }
    ifs.close();
    ofs.close();
}

void Graph::LoadUpdateStream(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);


    std::string type;
    while (ifs >> type)
    {
        if (type == "v" || type == "-v")
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            updates_.emplace('v', type == "v", vertex_id, 0u, label,0.0);
        }
        else
        {
            uint from_id, to_id, label;
            float weight;
            label=0;
            ifs >> from_id >> to_id >>label>>weight;
            updates_.emplace('e', type == "e", from_id, to_id, label,weight);
        }
    }
    ifs.close();
}

void Graph::PrintMetaData() const
{
    std::cout << "# vertices = " << NumVertices() <<
        "\n# edges = " << NumEdges() << std::endl;
}
void Graph::InitLabelIndex() {
    uint vNum=this->NumVertices();
    uint vlabelNum=this->NumVLabels();
    uint elabelNum= this->NumELabels();
    this->labelIndex.resize(vNum, nullptr);
    for(uint k=0;k<vNum;k++){
        std::vector<Neighbor>neighbors=vNeighbors[k];
        int* labelDistribution=new int[vlabelNum+elabelNum]();
        for(uint i=0;i<neighbors.size();i++){
            labelDistribution[neighbors[i].getVertexLabel()]++;
            labelDistribution[vlabelNum+neighbors[i].GetEdgelabel()]++;
        }

       this->labelIndex[k]=labelDistribution;
    }
}
void Graph::InitMatchOrderType(const std::vector<std::vector<uint> > &order_vs_,const std::vector<std::vector<std::vector<uint>>>&rightNeighbor) {
    //标记所有的查询顶点为自由匹配点和/或者孤立顶点
    for(uint i=0;i<this->NumEdges();i++){
        std::vector<vertexType> currentMatchOrderTypes;
        std::vector<std::vector<int>>currentMatchOrderLDvertex;
        currentMatchOrderLDvertex.resize(this->NumVertices(),{});
        for(int j=0;j<order_vs_[i].size();j++){
            int id=order_vs_[i][j];
            vertexType type;
            if(rightNeighbor[i][id].size()==0)
            {
                type=isolatedVertex;
            }
            else{
                type=freeVertex;
            }
           /* if(j==0||j==1){
                type=freeVertex;
            }
            else{
                for(int k=j+1;k<order_vs_[i].size();k++){
                    const auto &forwardneighbor=forwardNeighbors[i][k];
                    for(int m=0;m<forwardneighbor.size();m++){
                        if(forwardneighbor[m].GetVetexIndex()==j){
                                type=freeVertex;
                            break;
                        }
                    }
                    if(type==freeVertex)
                        break;
                }
            }*/
            currentMatchOrderTypes.emplace_back(type);
        }
        this->matchVertexTypes.push_back(currentMatchOrderTypes);
//        this->LDRecord.push_back(currentMatchOrderLDvertex);
    }


    //收集孤立顶点集合
    for(int i = 0; i < this->NumEdges(); i++){
        const auto & vertexTypes = this->matchVertexTypes[i];
        std::vector<int>isolateRecord;
        for(int k = 0; k < vertexTypes.size(); k++){
            if(vertexTypes[k] == isolatedVertex){
                isolateRecord.push_back(k);
            }
        }
        this->isolatedRecord.push_back(isolateRecord);
    }
    std::cout<<"Initial MatchOrderType end"<<std::endl;
}
void Graph::UpdateLabelIndex(uint v1, uint v2, uint label, uint flag) {
    uint v1Label = this->GetVertexLabel(v1);
    uint v2Label = this->GetVertexLabel(v2);
    uint Label = this->NumVLabels() + label;
    if(flag == 1){
        this->labelIndex[v1][v2Label]++;
        this->labelIndex[v1][Label]++;
        this->labelIndex[v2][v1Label]++;
        this->labelIndex[v2][Label]++;
    }
    else{
        this->labelIndex[v1][v2Label]--;
        this->labelIndex[v1][Label]--;
        this->labelIndex[v2][v1Label]--;
        this->labelIndex[v2][Label]--;
    }
}
const vertexType Graph::GetVertexType(uint order_index, uint depth) {
    return this->matchVertexTypes[order_index][depth];
}
std::vector<uint>Graph::GetIsolateVertexBeforeDepth(uint order_index, uint depth) {
    std::vector<uint> result;
    auto &isolatedVertexs=isolatedRecord[order_index];
    for(int i=0;i<isolatedVertexs.size();i++){
        if(isolatedVertexs[i]<depth){
            result.emplace_back(isolatedVertexs[i]);
        }
    }
    return result;
}
void Graph::setBatchVertexType(uint order_index, const std::vector<uint> &vertexs, vertexType type) {
    for(int i = 0; i < vertexs.size(); i++){
        this->matchVertexTypes[order_index][vertexs[i]] = type;
    }
}
bool Graph::isNeighbor(uint v1, uint v2) {
    const std::vector<uint> *nbrs;
    uint other;
    //从邻居少的节点来找边标签
    if (GetDegree(v1) < GetDegree(v2))
    {
        nbrs = &GetNeighbors(v1);
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        other = v1;
    }

    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)
    {
        mid = (start + end) / 2;
        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            return true;
        }
    }
    return false;
}