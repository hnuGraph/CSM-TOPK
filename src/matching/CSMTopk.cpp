#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <set>
#include <map>
#include "../utils/types.h"
#include "../utils/globals.h"
#include "../utils/utils.h"
#include "../graph/graph.h"
#include "CSMTopk.h"

struct pairCompare {
    bool operator()(const std::pair<float, int>& p1, const std::pair<float, int>& p2) {
        if (p1.first == p2.first) {
            return p1.second > p2.second;
        }
        return p1.first < p2.first;
    }
};


bool ForwardNeighborcmp(ForwardNeighbor*f1,ForwardNeighbor*f2){
    return (*f1)>(*f2);
}

CSMTopk::CSMTopk(Graph& query_graph, Graph& data_graph,Subgraph& global_subgraph,
                     uint max_num_results,
                     bool print_prep,
                     bool print_enum,
                     bool homo)
        : matching(query_graph, data_graph, global_subgraph,max_num_results,
                   print_prep, print_enum, homo)
        , order_vs_(query_.NumEdges())
        , order_csrs_(query_.NumEdges())
        , order_offs_(query_.NumEdges())
        ,order_vertex_index(query_.NumEdges())
        ,topKSet(0)
        ,suffixMax(query_graph.NumVertices(),0)
        ,isolatedMax(query_graph.NumVertices(),-1)
        ,rightNeighbor(query_graph.NumEdges())
        ,matchCandidate(query_graph.NumVertices())
        ,match(query_graph.NumVertices())
        ,labelToQueryVertex(query_graph.NumVLabels())
        ,globalVkMatchUk(data_graph.NumVertices())
        ,globalStarIndex(query_.NumEdges())
        ,queryVertexIndexInlabel(query_.NumVertices())
        ,LocalStarIndex(query_.NumVertices())
        ,matchLeftNeighborSum(query_.NumEdges())
        ,leftNeighborIdSum(query_.NumEdges())
        ,matchVetexLeftNeighbor(query_.NumVertices())
        , matchVetexSumweight(query_.NumVertices())
{
//    globalStarIndex.resize(query_.NumEdges());
    for (uint i = 0; i < query_.NumEdges(); ++i)
    {
        order_vs_[i].resize(query_.NumVertices());//节点个数
        order_csrs_[i].resize(query_.NumEdges() + 1);//边的个数+1，
        order_offs_[i].resize(query_.NumVertices(), 0);//节点个数，初始化为0
        globalStarIndex[i].resize(query_.NumVertices());
        order_vertex_index[i].resize(query_.NumVertices());
        rightNeighbor[i].resize(query_.NumVertices());
        matchLeftNeighborSum[i].resize(query_.NumVertices());
        leftNeighborIdSum[i].resize(query_.NumVertices());
    }
}
CSMTopk::~CSMTopk() noexcept {
    for(MatchRecord* item:topKSet){
        delete item;
        item= nullptr;
    }
    for(int i=0;i<query_.NumEdges();i++){
        for(StarGraph*s:globalStarIndex[i]){
            delete s;
        }
    }
}
/**
 * Preprocessing process
 */
void CSMTopk::Preprocessing()
{
    this->data_.InitLabelIndex();
    total_index_time.StartTimer();
    GenerateMatchingOrder();
    total_index_time.StopTimer();
    this->query_.InitLabelIndex();
    //create globalsubgraph;
    createGlobalSubgraph();

#ifdef LOG_TRACK
    stringstream _ss;
    for(int i=0;i<query_.NumVertices();i++){
        _ss<<i<<"candidate:"<<endl;
        for(auto m:globalsubgraph_.matchCandidate[i])
        {
            _ss<<m<<" ";
        }
        _ss<<endl;
    }
    Log::track1(_ss);
#endif
    this->query_.InitMatchOrderType(this->order_vs_,this->rightNeighbor);
    createLabelToQueryVertex();
    //createStarIndex with globalSubgraph
    CreateStarIndex();
    std::cout<<"Preprocess end"<<endl;
}
/**
 * update mwstar
 * @param match_index
 * @param caddidate_v
 * @param candidate_u
 * @param candidate_v_index
 */
void CSMTopk::updateStarIndex(uint match_index, uint caddidate_v, uint candidate_u,int candidate_v_index) {
    std::vector<int>&result=globalVkMatchUk[caddidate_v][match_index];
    int vertex_index=order_vertex_index[match_index][candidate_u];
    if(vertex_index==0)
    {
        return;
    }
    StarGraph* s=globalStarIndex[match_index][vertex_index];
    const std::vector<ForwardNeighbor*>&queryVetex=s->GetqueryVertex();
    std::vector<Neighbor>&vN= this->globalsubgraph_.vNeighbors[caddidate_v];
    int leftvN=0;
    int rightqV=0;
    int vNsize=vN.size();
    int qVsize=queryVetex.size();
    int flag=1;
    float sumweight=0;
   // std::vector<uint>MatchId;
    while(leftvN<vNsize&&rightqV<qVsize){
        if(vN[leftvN].GetelabelAndVertexLabel()<queryVetex[rightqV]->GetelabelAndVertexLabel())
        {
            flag=0;
            break;
        }

        while(vN[leftvN].GetelabelAndVertexLabel()>queryVetex[rightqV]->GetelabelAndVertexLabel()||vN[leftvN].getMatchQueryVertexId()!=queryVetex[rightqV]->GetVetexId()||vN[leftvN].getfromVertexId()!=candidate_u)
        {
            leftvN++;
            if(leftvN>=vN.size())
            {
                flag=0;
                break;
            }
        }
        if(!flag)
            break;
        if(vN[leftvN].GetelabelAndVertexLabel()==queryVetex[rightqV]->GetelabelAndVertexLabel()){
           // MatchId.emplace_back(vN[leftvN].getVertexId());
            float edgeweight=vN[leftvN].GetEdgeWeight();
            rightqV++;
            leftvN++;
            sumweight+=edgeweight;
        }
    }
    if(!flag)
    {
        return;
    }
    else if(rightqV==qVsize){
        //globalVkMatchUk update
        result[candidate_v_index]=sumweight;
        // globalStarIndex update
        if(s->getStarMaxWeight()==queryVetex.size()*mw||s->getStarMaxWeight()<sumweight)
        {
            s->setStarMaxWeight(sumweight);
            s->setMatchDataVertexId(caddidate_v);
        }
    }
}
float CSMTopk::GetBackWeight(uint order_index,uint depth) {
    float sum=0;
    uint n=query_.NumVertices();
    std::vector<uint>& matchOrder= this->order_vs_[order_index];
    for(int i=depth;i<n;i++){
        sum+=LocalStarIndex[i];
    }
    return sum;
}
/**
 * Initial mwstar
 */
void CSMTopk::CreateStarIndex() {
    int m=query_.NumEdges();
    int n=query_.NumVertices();
    for(int i=0;i<n;i++){
        const std::vector<uint>&candidate=globalsubgraph_.matchCandidate[i];
        const std::vector<uint>&candidate_u=labelToQueryVertex[query_.GetVertexLabel(i)];
        for(uint v:candidate){
            globalVkMatchUk[v].resize(m);
            for(int j=0;j<m;j++){
                globalVkMatchUk[v][j].resize(candidate_u.size());
                if(i==order_vs_[j][0])
                    continue;
                int candidate_index=queryVertexIndexInlabel[i];
                updateStarIndex(j,v,i,candidate_index);
            }

        }
        }
    }

 vector<int> CSMTopk::EdgeisInMatchOrder(uint v1, uint v2, uint v1label, uint v2label,uint velabel) {
    vector<int>result;
    for(int i=0;i<order_vs_.size();i++){
        uint u1=order_vs_[i][0];
        uint u2=order_vs_[i][1];
        uint u1label=query_.GetVertexLabel(u1);
        uint u2label=query_.GetVertexLabel(u2);
        uint qlabel=std::get<2>(query_.GetEdgeLabel(u1,u2));
        if((v1label==u1label&&v2label==u2label&&velabel==qlabel)||(v1label==u2label&&v2label==u1label&&velabel==qlabel))
        {
            result.emplace_back(i);
        }
    }
    return result;
}


void CSMTopk::GenerateMatchingOrder()
{
    // generate the initial matching order, order_*s_[0]
    std::vector<bool> visited(query_.NumVertices(), false);
    uint max_degree = 0u;
    //首先找到的是度最大的节点
    for (size_t i = 0; i < query_.NumVertices(); i++)
    {
        if (query_.GetDegree(i) > max_degree)
        {
            max_degree = query_.GetDegree(i);
            order_vs_[0][0] = i;
            order_vertex_index[0][i]=0;
        }
    }
    visited[order_vs_[0][0]] = true;

    // loop over all remaining positions of the order
    for (uint i = 1; i < query_.NumVertices(); ++i)
    {
        uint max_adjacent = 0;
        uint max_adjacent_u = NOT_EXIST;
        //Find the vertex that is not in the sequence but has the highest number of neighbors in the sequence
        for (size_t j = 0; j < query_.NumVertices(); j++)
        {
            uint cur_adjacent = 0u;
            if (visited[j]) continue;

            auto& q_nbrs = query_.GetNeighbors(j);
            for (auto& other : q_nbrs)
                if (visited[other])
                    cur_adjacent++;

            if (!cur_adjacent) continue;
            if (
                    max_adjacent_u == NOT_EXIST ||
                    (cur_adjacent == max_adjacent &&
                     query_.GetDegree(j) > query_.GetDegree(max_adjacent_u)) ||
                    cur_adjacent > max_adjacent
                    ) {
                max_adjacent = cur_adjacent;
                max_adjacent_u = j;
            }
        }
        order_vs_[0][i] = max_adjacent_u;
        order_vertex_index[0][max_adjacent_u]=i;
        visited[max_adjacent_u] = true;
        order_offs_[0][i] = order_offs_[0][i - 1];
        auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
        StarGraph*s=new StarGraph();
        for (auto &other: q_nbrs)
        {
            if (visited[other])
            {
                uint qlabel= std::get<2>(this->query_.GetEdgeLabel(max_adjacent_u,other));
                //  std::cout<<"globalStarIndex[0]["<<i<<"] other:"<<other<<" other label"<< this->query_.GetVertexLabel(other)<<endl;
                ForwardNeighbor* forwardNeighbor=new ForwardNeighbor(other, this->query_.GetVertexLabel(other),qlabel);
                s->AddForwardNeighbor(forwardNeighbor);
                order_csrs_[0][order_offs_[0][i]++] = other;
            }
        }
        s->InitalmaxWeight();
        globalStarIndex[0][i]=s;

    }

    // generate other incremental matching orders
    for (uint i = 1; i < query_.NumEdges(); ++i)
    {
        std::vector<bool> visited(query_.NumVertices(), false);

        // get the first edge
        std::vector<uint>::iterator it = std::lower_bound(
                order_offs_[0].begin(), order_offs_[0].end(), i + 1
        );
        uint tmp= *(order_vs_[0].begin() + std::distance(order_offs_[0].begin(), it));
        order_vs_[i][0]=tmp;
        order_vertex_index[i][tmp]=0;
        order_vs_[i][1] = order_csrs_[0][i];
        order_vertex_index[i][order_csrs_[0][i]]=1;
        order_csrs_[i][0] = order_vs_[i][0];
        StarGraph*s=new StarGraph();
        uint qlabel=std::get<2>(this->query_.GetEdgeLabel(order_vs_[i][0],order_vs_[i][1]));

        ForwardNeighbor *forwardNeighbor=new ForwardNeighbor(order_vs_[i][0], this->query_.GetVertexLabel(order_vs_[i][0]),qlabel);
        s->AddForwardNeighbor(forwardNeighbor);
        s->InitalmaxWeight();
        globalStarIndex[i][1]=(s);

        visited[order_vs_[i][0]] = true;
        visited[order_vs_[i][1]] = true;

        order_offs_[i][2] = order_offs_[i][1] = 1;
        for (uint j = 2; j < query_.NumVertices(); ++j)
        {
            uint max_adjacent = 0;
            uint max_adjacent_u = NOT_EXIST;
            for (size_t k = 0; k < query_.NumVertices(); k++)
            {
                uint cur_adjacent = 0u;
                if (visited[k]) continue;

                auto& q_nbrs = query_.GetNeighbors(k);
                for (auto& other : q_nbrs)
                    if (visited[other])
                        cur_adjacent++;

                if (!cur_adjacent) continue;
                if (
                        max_adjacent_u == NOT_EXIST ||
                        (cur_adjacent == max_adjacent &&
                         query_.GetDegree(k) > query_.GetDegree(max_adjacent_u)) ||
                        cur_adjacent > max_adjacent
                        ) {
                    max_adjacent = cur_adjacent;
                    max_adjacent_u = k;
                }
            }
            order_vs_[i][j] = max_adjacent_u;
            order_vertex_index[i][max_adjacent_u]=j;
            visited[max_adjacent_u] = true;

            order_offs_[i][j] = order_offs_[i][j - 1];
            StarGraph*s=new StarGraph();
            auto& q_nbrs = query_.GetNeighbors(max_adjacent_u);
            for (auto &other: q_nbrs)
            {
                if (visited[other])
                {
                    // std::cout<<"globalStarIndex["<<i<<"]"<<"["<<j<<"] "<<"other:"<<other<<" other label"<< this->query_.GetVertexLabel(other)<<endl;
                    order_csrs_[i][order_offs_[i][j]++] = other;
                    qlabel=std::get<2>(this->query_.GetEdgeLabel(max_adjacent_u,other));
                    ForwardNeighbor* forwardNeighbor=new ForwardNeighbor(other, this->query_.GetVertexLabel(other),qlabel);
                    s->AddForwardNeighbor(forwardNeighbor);
                }
            }
            s->InitalmaxWeight();
            globalStarIndex[i][j]=(s);
        }
    }
    //Sort global StarIndex according to <el, vl>for all matching orders
    for(int i=0;i<query_.NumEdges();i++){
        for(int j=1;j<query_.NumVertices();j++) {
            StarGraph *s = globalStarIndex[i][j];
            std::vector<ForwardNeighbor*> &globalIndex = s->GetqueryVertex();
            std::sort(globalIndex.begin(), globalIndex.end(), ForwardNeighborcmp);
        }
    }


    size_t sumId=0;
    //Create a right neighbor array
    if (print_preprocessing_results_)
    {
        std::cout << "matching order: " << std::endl;
        std::cout << "-vertex(backward neighbors)-\n";
        for (uint i = 0; i < query_.NumEdges(); ++i)
        {
            std::cout << "#" << i << ": ";
            for (uint j = 0; j < query_.NumVertices(); ++j)
            {
                sumId=0;
                std::cout << order_vs_[i][j];
                if (j == 0)
                {
                    std::cout << "-";
                    continue;
                }
                std::vector<ForwardNeighbor>currentQueryNeighbors;
                matchLeftNeighborSum[i][j]=order_offs_[i][j]-order_offs_[i][j-1];
                for (uint k = order_offs_[i][j - 1]; k < order_offs_[i][j]; k++)
                {
                    uint toVertexId=order_csrs_[i][k];
                    uint toVertexIndex=order_vertex_index[i][toVertexId];
                    uint toVertexLabel=query_.GetVertexLabel(order_csrs_[i][k]);
                    uint edgelabel=std::get<2>(query_.GetEdgeLabel(order_vs_[i][j],toVertexId));
                    ForwardNeighbor f(toVertexIndex,toVertexId,toVertexLabel,edgelabel);
                    currentQueryNeighbors.push_back(f);
                    rightNeighbor[i][toVertexId].emplace_back(order_vs_[i][j]);
                    if (k == order_offs_[i][j - 1]) std::cout << "(";
                    std::cout << order_csrs_[i][k];
                    if (k != order_offs_[i][j] - 1) std::cout << ",";
                    else std::cout << ")";
                }
                leftNeighborIdSum[i][j]=sumId;
                if (j != query_.NumVertices() - 1)
                    std::cout << "-";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}
void CSMTopk::InitialTopK(const std::string &path) {

    if (!io::file_exists(path.c_str()))
    {
        std::fstream fp(path,std::ios::out);
        for(auto t:topKSet){
            fp << t->printMatchRecord();
        }
        fp.close();
    }


#ifdef RESULT_TRACK
    stringstream _ss1;
    _ss1<<"Initial Top k"<<std::endl;
    // std::sort(topKSet.begin(),topKSet.end(), matchRecordCmp);
    for(auto d:topKSet){
        if(d!=NULL){
            _ss1<<d->toString();
            Log::track2(_ss1);
            _ss1.clear();
            _ss1.str("");

        }
    }
#endif
}

void CSMTopk::updateTopK() {
#ifdef PRINT_DEBUG
    /*for(auto it=edgeFlags.begin();it!=edgeFlags.end();it++){
        std::cout<<"edgeFlags["<<it->first.first<<","<<it->first.second<<"] "<<edgeFlags[std::make_pair(it->first.first,it->first.second)]<<std::endl;
    }*/
#endif
    stringstream _ss;
    if(!isUpdateIntopkset){
        return;
    }
    if(print_result){
        std::cout<<"after insert "<<std::endl;
        for(auto d:topKSet){
            std::cout<<d->toString();
        }
    }
#ifdef RESULT_TRACK
    _ss<<"after insert "<<std::endl;
    for(auto d:topKSet){
        _ss<<d->toString();
        Log::track2(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif
}

void CSMTopk::searchMatches(int depth,uint matchorderindex, searchType flag)  {
    //1.找到前向邻居，找到候选解
    std::vector<uint>& matchOrder= this->order_vs_[matchorderindex];
    uint queryVertex=matchOrder[depth];
    uint queryVertexLabel=this->query_.GetVertexLabel(queryVertex);
    std::vector<SingleCandidate>& singleVertexCandidate=this->matchCandidate[depth];
    std::vector<SingleCandidate>copySingleVertexCandidate=this->matchCandidate[depth];
    getIntersetSingleCandidate(singleVertexCandidate,matchorderindex,depth);
    if(singleVertexCandidate.size()==0){
        this->matchCandidate[depth]=copySingleVertexCandidate;
        return;
    }

    total_densityFilter_time.StartTimer();
    densityFilter(matchorderindex,depth,singleVertexCandidate);
    total_densityFilter_time.StopTimer();
    if(singleVertexCandidate.size()==0)
    {
        this->matchCandidate[depth]=copySingleVertexCandidate;
        return;
    }
    if(isInsert)
        IsearchSpace+=singleVertexCandidate.size();
    else
        DsearchSpace+=singleVertexCandidate.size();
    //Print_Time2("densityFilter ",start);
    //顺序扩展
    if(depth==query_.NumVertices()-1){
        //add matchresult;
        std::sort(singleVertexCandidate.begin(),singleVertexCandidate.end());
        for(const SingleCandidate &single:singleVertexCandidate){
            uint dataV=single.getVertexId();
            if(visited_[dataV])
                continue;
            float sumWeight= this->match[depth-1].getSumWeight();
            this->match[depth].setVertexId(dataV);
            sumWeight+=single.getSumWeight();
            this->match[depth].setSumWeight(sumWeight);
            int n=query_.NumVertices();
            std::vector<uint>m(n);
            for(int i=0;i<match.size();i++){
                m[order_vs_[matchorderindex][i]]=match[i].getVertexId();
            }
            float density=sumWeight/ (sqrt(n)*(n-1));

            MatchRecord *record=new MatchRecord(density,m);
            int matchResult= addMatchRecords(record);
            allMatchFind++;
            if(matchResult==1){
                if(flag==positive)
                {
                    num_positive_results_++;
                    numAddTopk++;
                }
            }
            else if(matchResult==3){
                this->matchCandidate[depth]=copySingleVertexCandidate;
                this->match[depth].clearSingleCandidate();
                return;
            }
        }
        //clear candidate;
        this->matchCandidate[depth]=copySingleVertexCandidate;
        this->match[depth].clearSingleCandidate();
        return;
    }
    else{
        for(int i=0;i<singleVertexCandidate.size();i++){
            uint dataV=singleVertexCandidate[i].getVertexId();
            if(visited_[dataV])
                continue;
            float weight=singleVertexCandidate[i].getSumWeight();
            //递归
            matchVertex(0,depth,dataV,weight);

            const std::vector<uint>&uk_neighbor=rightNeighbor[matchorderindex][queryVertex];
            std::vector<std::vector<SingleCandidate>>copyCandidate(uk_neighbor.size());
            std::vector<int>copyLocalStarIndex(query_.NumVertices());
            for(int i=0;i<uk_neighbor.size();i++){
                int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                copyCandidate[i]=matchCandidate[uk_neighbor_index];
                copyLocalStarIndex[uk_neighbor_index]=LocalStarIndex[uk_neighbor_index];
            }
            // std::cout<<"depth :"<<depth<<" data:"<<dataV<<endl;
            total_updaterightNeighborCandidate_time.StartTimer();
            bool isNull=updaterightNeighborCandidate(matchorderindex,queryVertex,0, false,dataV,uk_neighbor);
            total_updaterightNeighborCandidate_time.StopTimer();
            if(isNull)
            {
                for(int i=0;i<uk_neighbor.size();i++){
                    int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                    matchCandidate[uk_neighbor_index]= copyCandidate[i];
                    LocalStarIndex[uk_neighbor_index]=copyLocalStarIndex[uk_neighbor_index];
                }
                this->visited_[dataV]= false;
                continue;
            }
            //copy SingleCandidate
            //updateweight;
            searchMatches(depth+1,matchorderindex,flag);
            //返回candidate状态
            for(int i=0;i<uk_neighbor.size();i++){
                int uk_neighbor_index=order_vertex_index[matchorderindex][uk_neighbor[i]];
                matchCandidate[uk_neighbor_index]= copyCandidate[i];
                LocalStarIndex[uk_neighbor_index]=copyLocalStarIndex[uk_neighbor_index];
            }
            this->visited_[dataV]= false;
        }
        this->matchCandidate[depth]=copySingleVertexCandidate;
        this->match[depth].clearSingleCandidate();
    }
}



//flag==0 initial flag=1 update
void CSMTopk::FindMatches(uint flag,uint order_index, uint depth, std::vector<uint> m, size_t &num_results, float density_s) {
    if (reach_time_limit) return;
#ifdef PRINT_DEBUG
    if(density_s<0){
            std::cout<<"density <0 u:"<<order_vs_[order_index][depth]<<" depth: "<<depth<<"density: "<<density_s<<" match: ";
            for(int i=0;i<m.size();i++){
                std::cout<<m[i]<<" ";
            }
            std::cout<<std::endl;
        }
#endif
    if(flag==1) {
        float back_max_result = GetBackWeight(order_index, depth);
        uint n = query_.NumVertices();
        if (topKSet.size() == k) {
            float weight = topKSet.back()->getDensity();
            uint tmpdensity = density_s + back_max_result;
            if (tmpdensity / (sqrt(n) * (n - 1)) < weight)
                return;
        }
    }
    uint u = order_vs_[order_index][depth];
    uint u_min = NOT_EXIST;
    uint u_min_label = NOT_EXIST;
    uint u_min_size = UINT_MAX;


    // find u_min
    const auto &q_nbrs = query_.GetNeighbors(u);
    const auto &q_nbr_labels = query_.GetNeighborLabels(u);
    for (uint i = 0u; i < q_nbrs.size(); i++) {
        //q_nbrs=0,3  u_other=0
        const uint u_other = q_nbrs[i];
        const uint u_other_label = q_nbr_labels[i];

        if (m[u_other] == UNMATCHED) continue;

        const uint cur_can_size = data_.GetNeighbors(m[u_other]).size();
        if (cur_can_size < u_min_size) {
            u_min_size = cur_can_size;
            u_min = u_other;
            u_min_label = u_other_label;
        }
    }

    const auto &u_min_nbrs = data_.GetNeighbors(m[u_min]);
    const auto &u_min_nbr_labels = data_.GetNeighborLabels(m[u_min]);

    float tmp;
    bool candidate_empty = true;
    for (uint i = 0u; i < u_min_nbrs.size(); i++) {
        const uint v = u_min_nbrs[i];
        tmp = 0;
        // 1. check labels
        num_intermediate_results_before_index_check_++;
        if (
                data_.GetVertexLabel(v) != query_.GetVertexLabel(u) ||
                u_min_nbr_labels[i] != u_min_label
                )
            continue;
        num_intermediate_results_after_index_check_++;


        tmp += data_.GetEdgeWeight(m[u_min], v);


        // 2. check if joinable
        bool joinable = true;
        for (uint j = 0u; j < q_nbrs.size(); j++) {
            const uint u_other = q_nbrs[j];
            const uint u_other_labels = q_nbr_labels[j];

            if (m[u_other] == UNMATCHED || u_other == u_min) continue;

            auto it = std::lower_bound(data_.GetNeighbors(m[u_other]).begin(), data_.GetNeighbors(m[u_other]).end(),
                                       v);
            uint dis = std::distance(data_.GetNeighbors(m[u_other]).begin(), it);
            if (
                    it == data_.GetNeighbors(m[u_other]).end() ||
                    *it != v ||
                    data_.GetNeighborLabels(m[u_other])[dis] != u_other_labels
                    ) {
                joinable = false;
                break;
            }
            tmp += data_.GetEdgeWeight(m[u_other], v);
        }
        if (!joinable) continue;
        num_intermediate_results_after_joinability_check_++;

        candidate_empty = false;

        // 3. check if visited
        if (!homomorphism_ && visited_[v]) continue;
        num_intermediate_results_after_visit_check_++;

        // 4. add a vertex mapping

        m[u] = v;
        visited_[v] = true;
        density_s += tmp;

        if (depth == query_.NumVertices() - 1) {
            float lastds = density_s / (sqrt(m.size()) * (m.size() - 1));
            //sort(m.begin(),m.end());
            num_results++;
            MatchRecord *r = new MatchRecord(lastds, m);
            addMatchRecords(r);


            if (print_enumeration_results_) {
                for (auto j: m) {

                    std::cout << j << " ";
                }
            }

        } else {
            size_t num_results_before_recursion = num_results;
            FindMatches(flag,order_index, depth + 1, m, num_results, density_s);
            if (num_results == num_results_before_recursion) {
                num_intermediate_results_without_results_++;
            }

        }
        visited_[v] = false;
        m[u] = UNMATCHED;


#ifdef PRINT_DEBUG
        if (density_s < tmp) {
                std::cout << "u:" << order_vs_[order_index][depth] << " depth: " << depth << "density: " << density_s
                          << " tmp:" << tmp << " match: ";
                for (int i = 0; i < m.size(); i++) {
                    std::cout << m[i] << " ";
                }
                std::cout << std::endl;
            }
#endif
        density_s -= tmp;

//            tmin = copytmin;//回溯将tmin转为初始状态
        if (num_results >= max_num_results_) return;
        if (reach_time_limit) return;
    }
}


int CSMTopk:: addMatchRecords(MatchRecord *r) {
    int n=topKSet.size();
    if(n<k){
        for(int j=n-1;j>=0;j--){
            if((*topKSet[j])==(*r)){
                delete r;
                return 2;
            }
        }
        for(int j=n-1;j>=0;j--){
            if((*topKSet[j])>(*r)){
                topKSet.insert(topKSet.begin()+j+1,r);
                break;
            }
            if(j==0){
                topKSet.insert(topKSet.begin(),r);
            }
        }
        if(n==0){
            topKSet.insert(topKSet.begin(),r);
        }
        isUpdateIntopkset= true;
        return 1;
    }
    else{
        for(int j=n-1;j>=0;j--){
            if((*topKSet[j])==(*r)){
                delete r;
                return 2;
            }
        }
        MatchRecord* d = topKSet.back();
        if((*r)>(*d))
        {
            delete d;
            topKSet.pop_back();
            int m=topKSet.size();
            for(int j=m-1;j>=0;j--){
                if((*topKSet[j])>(*r)){
                    topKSet.insert(topKSet.begin()+j+1,r);
                    break;
                }
                if(j==0){
                    topKSet.insert(topKSet.begin(),r);
                }
            }
#ifdef PRINT_DEBUG
            stringstream _ss;
            _ss<<"insert record: ";
            _ss<<r->toString();
            _ss<<"after insert top k "<<std::endl;
            for(auto d:topKSet){
                _ss<<d->toString();
                Log::track1(_ss);
                _ss.clear();
                _ss.str("");
            }
#endif
            isUpdateIntopkset=true;
            return 1;
        }
        else{
            delete r;
            return 3;
        }
    }
}

void CSMTopk::InitialMatching(const std::string &path) {

    if (!io::file_exists(path.c_str()))
    {
        std::cout << "the file not exit " << path << std::endl;
        std::vector<uint> m(query_.NumVertices(), UNMATCHED);
        uint flag=0;
        float density_s=0;
        uint tmin=INT_MAX;
        uint order_index=0;
        uint depth=1;
        //#pragma omp parallel for num_threads(10) firstprivate(m) firstprivate(visited_) firstprivate(density_s) firstprivate(flag) firstprivate(tmin) firstprivate(order_index) firstprivate(depth)
        for (size_t i = 0; i < data_.NumVertices(); i++)
        {
            //std::cout<<"thread id"<<omp_get_thread_num<<endl;
            if (data_.GetVertexLabel(i) != NOT_EXIST)
            {
#ifdef PRINT_DEBUG
                stringstream _ss;
                _ss<<"vertex id:"<<i<<std::endl;
                Log::track1(_ss);
#endif
                if (query_.GetVertexLabel(order_vs_[0][0]) == data_.GetVertexLabel(i))
                {
                    m[order_vs_[0][0]] = i;
                    visited_[i] = true;

                    FindMatches(flag,order_index, depth, m, num_initial_results_,density_s);
                    visited_[i] = false;
                    m[order_vs_[0][0]] = UNMATCHED;
                }

            }

        }
    }
    else{
        std::ifstream ifs2(path);
        std::cout<<"load topk from file...."<<std::endl;
        char type;
        uint cnt=0;
        while (ifs2 >> type){
            if(cnt>=k)
                break;
            if(type == 't'){
                float density;
                uint tmp;
                std::vector<uint>m;
                ifs2>>density;
                for(int i=0;i<query_.NumVertices();i++){
                    ifs2>>tmp;
                    m.emplace_back(tmp);
                }
                MatchRecord* matchRecord=new MatchRecord(density,m);
                topKSet.push_back(matchRecord);
                //addMatchRecords(matchRecord);
            }
            cnt++;
        }
    }
}

//动态的加边操作
void CSMTopk::AddEdge(uint v1, uint v2, uint label, float weight) {
#ifdef LOG_TRACK
     stringstream _ss;

    for(auto n:globalsubgraph_.vNeighbors[63117]){
        if(n.getVertexId()==1709&&n.getMatchQueryVertexId()==0){
            _ss<<"find 1709"<<endl;
        }
    }
    Log::track1(_ss);
#endif
    total_update_globalIndex_time.StartTimer();
    int numAddTopk=0;
    allMatchFind=0;
    bool flag= true;
    data_.AddEdge(v1, v2, label, weight,  1);
    this->data_.UpdateLabelIndex(v1,v2,label,1);
    uint v1label=data_.GetVertexLabel(v1);
    uint v2label=data_.GetVertexLabel(v2);
    vector<int>match=EdgeisInMatchOrder(v1,v2,v1label,v2label,label);
    if(match.size()==0){
        total_update_globalIndex_time.StartTimer();
        return;
    }
    //update globalsubgraph and starIndex
    bool isInGlobalSubgraph=updateGlobalSubgraph(v1,v2,label,weight,match);;
    if(!isInGlobalSubgraph){
        total_update_globalIndex_time.StopTimer();
        return;
    }
    total_update_globalIndex_time.StopTimer();
#ifdef LOG_TRACK
    stringstream _ss1;
    for(auto n:globalsubgraph_.matchCandidate[1]){
        if(n==62501)
        {
            _ss<<"find 62501"<<endl;
        }
    }
#endif
    total_search_time.StartTimer();
    isUpdateIntopkset= false;
    for(auto m:match){
        uint u1=order_vs_[m][0];
        uint u2=order_vs_[m][1];
        uint u1label= this->query_.GetVertexLabel(u1);
        uint u2label=this->query_.GetVertexLabel(u2);
        uint v1label=this->data_.GetVertexLabel(v1);
        uint v2label=this->data_.GetVertexLabel(v2);
        float weight=this->data_.GetEdgeWeight(v1,v2);
        InitialLocalIndex(m);
        if(v1label!=v2label){
            if(v1label!=u1label)
            {
                swap(v1,v2);
            }
            //todo
           flag=SearchMatchesWithEdge(m,v1,v2,weight,u1,u2,positive);
            if(flag)
                continue;
        }
        else{
            for(int i = 0; i < 2; i++) {
              flag= SearchMatchesWithEdge(m,v1,v2,weight,u1,u2,positive);
                if(flag)
                {
                    std::swap(v1, v2);
                    continue;
                }
                std::swap(v1, v2);
            }
        }
    }
    total_search_time.StopTimer();
    //Print_Time2("SearchMatches ", start);
    END_ENUMERATION:
    total_print_time.StartTimer();
    sumAllMatchFind+=allMatchFind;
//    std::cout<<"num add top k:"<<numAddTopk<<endl;
//    std::cout<<"all match find:"<<allMatchFind<<endl;
    updateTopK();
    total_print_time.StopTimer();
//    Print_Time2("PrintTopk ", start);
}




void CSMTopk::deleteUpdateTopK() {
    if(print_result) {
        std::cout << "after delete" << std::endl;
        for (auto d: topKSet) {
            std::cout << d->toString();
        }
    }
#ifdef RESULT_TRACK
    stringstream _ss;
    _ss<<"after delete"<<std::endl;
    //std::sort(topKSet.begin(),topKSet.end(), matchRecordCmp);
    for(auto d:topKSet){
        _ss<<d->toString();
        Log::track2(_ss);
        _ss.clear();
        _ss.str("");
    }
#endif

}

//动态地减边操作
void CSMTopk::RemoveEdge(uint v1, uint v2,uint label) {
    //1.更新data、更新nlf标签
    total_delete_update_time.StartTimer();
    allMatchFind=0;
    uint v1label=data_.GetVertexLabel(v1);
    uint v2label=data_.GetVertexLabel(v2);
    float weight=data_.GetEdgeWeight(v1,v2);
    data_.RemoveEdge(0,v1, v2);
    data_.UpdateLabelIndex(v1,v2,label,0);
    bool flag;
    //2.更新subgraph中的点和边
    vector<int>match=EdgeisInMatchOrder(v1,v2,v1label,v2label,label);
    //3.删除边并且更新索引
    deleteGlobalSubgraph(v1,v2,label,weight,match);
    total_delete_update_time.StopTimer();
#ifdef LOG_TRACK
 /*   stringstream _ss;
    if(globalStarIndex[0][7]->getStarMaxWeight()==1){
        _ss<<"max=1"<<endl;
    }
    Log::track1(_ss);*/

#endif
    //4.判断topk中是否包含删除边
    total_delete_time.StartTimer();
     flag=deleteMatchRecordWithEdge(v1,v1label,v2,v2label,label,match);
    if(!flag)
    {
        total_delete_time.StopTimer();
        return;
    }
    //5 从subgraph中的Edge进行重搜
    //从候选最少的节点，以及它的邻居候选次少的节点进行重搜
    int n=query_.NumVertices();
    uint minVertexSize=UINT_MAX;
    uint minVertex=0;
    for(int i=0;i<n;i++){
        int s=globalsubgraph_.matchCandidate[i].size();
        if(s<minVertexSize){
            minVertexSize=s;
            minVertex=i;
        }
    }
    const std::vector<uint>&minNeighbor=query_.GetNeighbors(minVertex);
    uint minVertexNeighorSize=UINT_MAX;
    uint minVertexNeighor=0;
    for(uint u:minNeighbor){
        int us=globalsubgraph_.matchCandidate[u].size();
        if(us<minVertexNeighorSize){
            minVertexNeighorSize=us;
            minVertexNeighor=u;
        }
    }
    int m=query_.NumEdges();
    int matchIndex=0;
    for(int i=0;i<m;i++){
        if((order_vs_[i][0]==minVertex&&order_vs_[i][1]==minVertexNeighor)||((order_vs_[i][1]==minVertex&&order_vs_[i][0]==minVertexNeighor))){
            matchIndex=i;
            break;
        }
    }
    const std::vector<uint>&minVertexCandidate=globalsubgraph_.matchCandidate[minVertex];
    for(uint mv:minVertexCandidate){
        const std::vector<Neighbor>&neighbors=globalsubgraph_.GetVNeighbors(mv);
        for(const Neighbor&n:neighbors){
            if(n.getMatchQueryVertexId()==minVertexNeighor&&n.getfromVertexId()==minVertex){
                //search
                uint u1=order_vs_[matchIndex][0];
                uint u2=order_vs_[matchIndex][1];
                uint u1label= this->query_.GetVertexLabel(u1);
                uint u2label=this->query_.GetVertexLabel(u2);
                uint v1=mv;
                uint v2=n.getVertexId();
                uint v1label=this->data_.GetVertexLabel(mv);
                uint v2label=n.getVertexLabel();
                float weight=n.GetEdgeWeight();
                InitialLocalIndex(matchIndex);

                if(v1label!=v2label){
                    if(v1label!=u1label)
                    {
                        swap(v1,v2);
                    }
                    //todo
                    SearchMatchesWithEdge(matchIndex,v1,v2,weight,u1,u2,negative);

                }
                else{
                    for(int i = 0; i < 2; i++) {
                        SearchMatchesWithEdge(matchIndex,v1,v2,weight,u1,u2,negative);
                        std::swap(v1, v2);
                    }
                }
            }
        }
    }
    total_delete_time.StopTimer();
    sumDeleteallMatchFind+=allMatchFind;
//    std::cout<<"delete research matches:"<<allMatchFind<<endl;
//    std::cout<<"v1:"<<v1<<" v2:"<<v2<<endl;
    deleteUpdateTopK();
#ifdef PRINT_INDEX
    stringstream _ss;
    _ss<<(10000-curupdatenum)+1<<" "<<allMatchFind<<endl;
    Log::track3(_ss);
#endif
}

void CSMTopk::AddVertex(uint id, uint label) {
    data_.AddVertex(id, label);
    visited_.resize(id + 1, false);
}

void CSMTopk::RemoveVertex(uint id) {
    data_.RemoveVertex(id);
}

void CSMTopk::GetMemoryCost(size_t &num_edges, size_t &num_vertices) {
    num_edges = 0ul;
    num_vertices = 0ul;
}
bool CSMTopk::LabelFilter(uint data_v, uint query_v) {
    uint dataVlabelNum= this->data_.NumVLabels();
    uint dataElabelNum=this->data_.NumELabels();
    uint queryVlabelNum=this->query_.NumVLabels();
    uint queryElabelNum=this->query_.NumELabels();
    const auto & dataLabelIndex= this->data_.labelIndex[data_v];
    const auto & queryLabelIndex=this->query_.labelIndex[query_v];

    for(int i=0;i<queryVlabelNum;i++){
        if(i<dataVlabelNum){
            if(dataLabelIndex[i]<queryLabelIndex[i])
            {
                return false;
            }
        }
        else{
            if(queryLabelIndex[i]>0){
                return false;
            }
        }
    }
    for(int i=0;i<queryElabelNum;i++){
        if(i<dataElabelNum){
            if(dataLabelIndex[dataVlabelNum+i]<queryLabelIndex[queryVlabelNum+i]){
                return false;
            }
        }
        else{
            if(queryLabelIndex[queryVlabelNum+i]>0)
                return false;
        }
    }
    return true;
}




void CSMTopk::matchVertex(bool isFirstEdge,uint depth,uint data_v,float w) {
    if(isFirstEdge){
        this->match[depth].setVertexId(data_v);
        this->match[depth].setSumWeight(w);
        this->matchCandidate[depth].emplace_back(SingleCandidate(data_v,w));
        this->visited_[data_v]= true;
    }
    else{
        int preindex=0;
        for(int i=depth-1;i>=0;i--){
            if(match[i].getVertexId()!=-1){
                preindex=i;
                break;
            }
        }
        float weight=match[preindex].getSumWeight()+w;
        this->match[depth].setVertexId(data_v);
        this->match[depth].setSumWeight(weight);
        this->visited_[data_v]= true;
    }

}

void CSMTopk::popVertex(uint depth,uint data_v) {
    this->match[depth].clearSingleCandidate();
    this->matchCandidate[depth].clear();

    this->visited_[data_v]=false;
}
void CSMTopk::popVertex(uint data_v, uint matchorderindex, uint depth, const std::vector<uint>&uk_neighbor) {
    this->match[depth].clearSingleCandidate();
    const int n=uk_neighbor.size();
    for(int u_id:uk_neighbor){
        int query_order_index=order_vertex_index[matchorderindex][u_id];
        matchCandidate[query_order_index].clear();
    }
    this->matchCandidate[depth].clear();

    this->visited_[data_v]=false;
}

void CSMTopk::densityFilter(uint matchorder_index,uint depth,std::vector<SingleCandidate>&singleVertexCandidate) {
    uint n=query_.NumVertices();
    if(topKSet.size()<k){
        auto iter=singleVertexCandidate.begin();
        while(iter!=singleVertexCandidate.end()) {
            float md = (*iter).getSumWeight();
            iter++;
        }
        return;
    }
    float kw=topKSet.back()->getDensity();
    float sumWeight=0;
    //2. singleVertex prune
    bool flag= false;
    // extends
    sumWeight+=this->match[depth-1].getSumWeight();


    float backWeight=GetBackWeight(matchorder_index,depth+1);
    sumWeight+=backWeight;
    int cnt=0;

    auto iter=singleVertexCandidate.begin();
    while(iter!=singleVertexCandidate.end()){
        float md=(*iter).getSumWeight();
        float tmpweight=(sumWeight+md)/(sqrt(n)*(n-1));
        if(tmpweight<kw){
            (*iter).setVertexId(-1);
            cnt++;
        }
        iter++;
    }
    if(cnt==0)
        return;
    int newLen=singleVertexCandidate.size()-cnt;
    int svcLen=singleVertexCandidate.size();
    if(newLen==0)
    {
        singleVertexCandidate.resize(0);
        return;
    }
    else{
        int i=0;
        int j=0;
        while(j<svcLen){
            if(singleVertexCandidate[j].getVertexId()!=-1){
                singleVertexCandidate[i]=singleVertexCandidate[j];
                i++;
            }
            j++;

        }
        singleVertexCandidate.resize(newLen);
    }
}

void CSMTopk::createLabelToQueryVertex() {
    for(int i=0;i<query_.NumVertices();i++){
        uint label=query_.GetVertexLabel(i);
        labelToQueryVertex[label].emplace_back(i);
        queryVertexIndexInlabel[i]=labelToQueryVertex[label].size()-1;
    }
}

bool CSMTopk::updaterightNeighborCandidate(int matchorderindex,uint uk,uint uk_neigh,bool isFirstEdge, uint vk,const std::vector<uint>&uk_neighbor) {
    const  std::vector<Neighbor>&vN= this->globalsubgraph_.vNeighbors[vk];
#ifdef LOG_TRACK
#endif
    const int n=uk_neighbor.size();
    //1.find candidates of all right neighbors
    for(int i=0;i<n;i++) {
        uint query_id = uk_neighbor[i];
        if(isFirstEdge){
            if(query_id==uk_neigh)
            {
                isFirstEdge= false;
                continue;
            }
        }
        uint query_vertex_label = query_.GetVertexLabel(query_id);
        int query_order_index=order_vertex_index[matchorderindex][query_id];
        uint query_elabel=std::get<2>(query_.GetEdgeLabel(uk,query_id));
        StarGraph*s=globalStarIndex[matchorderindex][query_order_index];
        bool isFirstVertex= true;
        bool isCandidateFirstNull= true;
        float maxweight=0;
        float curWeight=0;
        if(matchCandidate[query_order_index].size()!=0)
            isCandidateFirstNull= false;
        //for vk's neighbor
        for (const auto &neighbor: vN) {
            uint neighbor_id = neighbor.getVertexId();
            const uint neighbor_match_id=neighbor.getMatchQueryVertexId();
            const uint neighbor_from_id=neighbor.getfromVertexId();
            if(visited_[neighbor_id]||neighbor_match_id!=query_id||neighbor_from_id!=uk)
                continue;
            uint  v_elabel=neighbor.GetEdgelabel();
            if (neighbor.getVertexLabel() == query_vertex_label&&query_elabel==v_elabel) {
                maxweight = globalVkMatchUk[neighbor_id][matchorderindex][queryVertexIndexInlabel[query_id]];
                //update LocalStarIndex
                //add candidate
                curWeight=neighbor.GetEdgeWeight();
                if (isCandidateFirstNull) {
                    if(isFirstVertex){
                        isFirstVertex= false;
                        LocalStarIndex[query_order_index]=maxweight;
                    }
                    else{
                        if(LocalStarIndex[query_order_index]<maxweight){
                            LocalStarIndex[query_order_index]=maxweight;
                        }
                    }
                    matchCandidate[query_order_index].emplace_back( neighbor_id,curWeight);
                }
                else{
                    std::vector<SingleCandidate>&neigh_candidate=matchCandidate[query_order_index];
                    for(SingleCandidate&s:neigh_candidate){
                        if(s.getVertexId()==neighbor_id){
                            if(isFirstVertex){
                                isFirstVertex= false;
                                LocalStarIndex[query_order_index]=maxweight;
                            }
                            else{
                                if(LocalStarIndex[query_order_index]<maxweight){
                                    LocalStarIndex[query_order_index]=maxweight;
                                }
                            }
                            s.addFlag();
                            s.addSumWeight(curWeight);
                            break;
                        }
                    }
                }
            }
            if(isInsert)
                IdeterminCandite++;
            else
                DdeterminCandite++;
        }
        if(matchCandidate[query_order_index].size()==0) {
            //isFirst recover
            for (int i = 0; i < n; i++) {
                uint query_id = uk_neighbor[i];
                int query_order_index = order_vertex_index[matchorderindex][query_id];
                matchCandidate[query_order_index].clear();
            }
            return true;
        }
    }
    return false;
}
void CSMTopk::InitialLocalIndex(int matchorderindex) {
    const std::vector<StarGraph*> & gs=globalStarIndex[matchorderindex];
    int n=gs.size();
    for(int i=1;i<n;i++){
        LocalStarIndex[i]=gs[i]->getStarMaxWeight();
    }
}
void CSMTopk::getIntersetSingleCandidate(std::vector<SingleCandidate> &singleVertexCandidate,int matchorderindex,int depth) {
    int i=0;
    int j=0;
    int csize=singleVertexCandidate.size();
    int len=0;
    int flag=matchLeftNeighborSum[matchorderindex][depth];
    if(flag==1)
        return;
    while(j<csize){
        if(singleVertexCandidate[j].getFlag()==flag){
            singleVertexCandidate[i]=singleVertexCandidate[j];
            i++;
            len++;
        }
        j++;

    }
    singleVertexCandidate.resize(len);
}
void CSMTopk::PrintAverageTime(int len) {
     int ilen=10000-len;
    std::cout <<"average query graph degree:"<< std::fixed << std::setprecision(2)<<query_.NumEdges()*2.0/query_.NumVertices()<<endl;
    std::cout<<"average data graph degree:"<<std::fixed << std::setprecision(2)<<data_.NumEdges()*2.0/data_.NumVertices()<<endl;
    std::cout << "average serach time: " << std::fixed << std::setprecision(2)
              << total_search_time.GetTimer() * 1.0 / ilen << " microseconds" << endl;
    std::cout << "average update global index time: " << std::fixed << std::setprecision(2)
              << total_update_globalIndex_time.GetTimer() * 1.0 / ilen << " microseconds" << endl;
    std::cout << "average update time " << std::fixed << std::setprecision(2)
              << (total_update_globalIndex_time.GetTimer() * 1.0 / ilen + total_search_time.GetTimer() * 1.0 / ilen)
              << " microseconds" << endl;
    std::cout << "average delete search time:" << std::fixed << std::setprecision(2)
              << total_delete_time.GetTimer() * 1.0 / len << " microseconds" << endl;
    std::cout << "average delete update global subgraph time:" << std::fixed << std::setprecision(2)
              << (total_delete_update_time.GetTimer() * 1.0) / len << " microseconds" << endl;
    std::cout << "average delete update time:" << std::fixed << std::setprecision(2)
              << (total_delete_time.GetTimer() * 1.0 / len + total_delete_update_time.GetTimer() * 1.0 / len)
              << " microseconds" << endl;
#ifdef COMPUTE_TRACK
stringstream _ss;
    _ss<< (total_update_globalIndex_time.GetTimer() * 1.0 / ilen + total_search_time.GetTimer() * 1.0 / ilen)<<","
       <<total_update_globalIndex_time.GetTimer() * 1.0 / ilen<<","
       <<total_search_time.GetTimer() * 1.0 / ilen<<","
       <<(total_delete_time.GetTimer() * 1.0 / len + total_delete_update_time.GetTimer() * 1.0 / len)<<","
       <<(total_delete_update_time.GetTimer() * 1.0) / len<<","
       <<total_delete_time.GetTimer() * 1.0 / len<<","
       <<sumAllMatchFind<<","
       <<sumDeleteallMatchFind<<","
       <<query_.NumEdges()*2.0/query_.NumVertices()<<","
       <<data_.NumEdges()*2.0/data_.NumVertices()<<","
       <<Itotal_densityfilter_time * 1.0 / ilen<<","
       <<(Itotal_updaterightNeighborCandidate_time* 1.0 / ilen)<<","
       <<total_densityFilter_time.GetTimer() * 1.0 / len<<","
       <<total_updaterightNeighborCandidate_time.GetTimer() * 1.0 / len<<","
       << IsearchSpace<<","
       <<DsearchSpace<<","
       <<IdeterminCandite<<","
       <<DdeterminCandite<<","
       <<space_cost<<endl;
    Log::track3(_ss);
#endif
}
void CSMTopk::createGlobalSubgraph() {
    //for q's vertices nlf check
    int n = data_.NumVertices();
    int m = query_.NumVertices();
    for (int i = 0; i < n; i++) {
        uint vlabel = data_.GetVertexLabel(i);
        for (int j = 0; j < m; j++) {
            uint qlabel = query_.GetVertexLabel(j);
            if (vlabel == qlabel) {
                if (LabelFilter(i, j)) {
                    globalsubgraph_.addQueryVertexCandidate(j, i);
                }
            }
        }
    }
    //addEdge
    for(int i=0;i< this->query_.NumVertices();i++){
        std::vector<uint>& irightneighbors=rightNeighbor[0][i];
        const std::vector<uint>&i_candidate=globalsubgraph_.matchCandidate[i];
        for(uint ic:i_candidate){
            const std::vector<uint>&i_nbrs=data_.GetNeighbors(ic);
            for(uint inbr:i_nbrs){
                for(uint j:irightneighbors){
                    const std::vector<uint>&j_candidate=globalsubgraph_.matchCandidate[j];
                    if(std::binary_search(j_candidate.begin(),j_candidate.end(), inbr)){
                        uint v1label=data_.GetVertexLabel(ic);
                        uint v2label=data_.GetVertexLabel(inbr);
                        uint edgelabel=std::get<2>(data_.GetEdgeLabel(ic,inbr));
                        float weight=data_.GetEdgeWeight(ic,inbr);
                        globalsubgraph_.AddEdge(i,j,ic,data_.GetVertexLabel(ic),inbr,data_.GetVertexLabel(inbr),edgelabel,weight);
                    }
                }
            }
        }
    }


}
bool CSMTopk::updateGlobalGraphHelp(int m, uint u1, uint u2, uint u1label, uint u2label, uint v1, uint v2, uint v1label, uint v2label,
                                    uint elabel, const std::vector<std::vector<uint>>&mcandidate, bool &flag) {
    bool isMatch= false;
    uint n=query_.NumEdges();
    bool isNew1= false;
    bool isNew2= false;
    bool isContain1= true;
    bool isContain2= true;
    if(v1label!=u1label)
    {
        swap(v1,v2);
    }
    if(!std::binary_search(mcandidate[u1].begin(),mcandidate[u1].end(), v1)){
        if(this->LabelFilter(v1,u1)){
            isNew1= true;
        }
        else{
            isContain1= false;
        }
    }
    if(!std::binary_search(mcandidate[u2].begin(),mcandidate[u2].end(), v2)){
        if(this->LabelFilter(v2,u2)){
            isNew2= true;
        }
        else{
            isContain2= false;
        }
    }
    if(isContain1&&isContain2){
        flag= true;
        isMatch= true;
        if(isNew1&&isNew2){
            updateglobalVertexStarIndex(u1,v1,u1label,elabel,n,mcandidate);
            updateglobalVertexStarIndex(u2,v2,u2label,elabel,n,mcandidate);
        }
        else if(isNew1&&!isNew2){
            updateglobalVertexStarIndex(u1,v1,u1label,elabel,n,mcandidate);
        }
        else if(!isNew1&&isNew2){
            updateglobalVertexStarIndex(u2,v2,u2label,elabel,n,mcandidate);
        }
        else{
            float w=data_.GetEdgeWeight(v1,v2);
            globalsubgraph_.AddEdge(u1,u2,v1,u1label,v2,u2label,elabel,w);
            int candidate_index= queryVertexIndexInlabel[u1];
            for(int j=0;j<query_.NumEdges();j++) {
                updateStarIndex(j,v1,u1,candidate_index);
            }
            int candidate_index2= queryVertexIndexInlabel[u2];
            for(int j=0;j<query_.NumEdges();j++) {
                updateStarIndex(j,v2,u2,candidate_index2);
            }
            numupdatestar+=2;
        }
    }
    else if(isContain1&&isNew1){
        updateglobalVertexStarIndex(u1,v1,u1label,elabel,n,mcandidate);
    }
    else if(isContain2&&isNew2){
        updateglobalVertexStarIndex(u2,v2,u2label,elabel,n,mcandidate);
    }
    return isMatch;
 }
bool CSMTopk::updateGlobalSubgraph(uint v1, uint v2, uint label, float weight,std::vector<int>&match) {
    bool flag= false;
        // 1.1Check if the nlf condition ?
        // 1.2true check if it is already in the candidate set.
        // 1.2.1 true, add an edge .
        // 1.2.2 false, update the candidate set
        uint n=query_.NumEdges();
        const auto &mcandidate=globalsubgraph_.matchCandidate;
        for(auto it=match.begin();it!=match.end();){
            bool isMatch= false;
            uint u1=order_vs_[*it][0];
            uint u2=order_vs_[*it][1];
            uint u1label=query_.GetVertexLabel(u1);
            uint u2label=query_.GetVertexLabel(u2);
            uint v1label=data_.GetVertexLabel(v1);
            uint v2label=data_.GetVertexLabel(v2);
            uint elabel=std::get<2>(query_.GetEdgeLabel(u1,u2));
            if(v1label!=v2label) {
                if(updateGlobalGraphHelp(*it, u1, u2, u1label, u2label, v1, v2, v1label, v2label, elabel,
                                      globalsubgraph_.matchCandidate, flag)){
                    isMatch= true;
                }
            }
            else{
                for(int i=0;i<2;i++){
                    if(updateGlobalGraphHelp(*it, u1, u2, u1label, u2label, v1, v2, v1label, v2label, elabel,
                                          globalsubgraph_.matchCandidate, flag)){
                        isMatch= true;
                    }
                    else{
                        if(isMatch){
                            continue;
                        }
                    }
                    std::swap(v1,v2);
                }
            }
           if(isMatch== false){
               it=match.erase(it);
           }
           else{
               ++it;
           }
        }
#ifdef LOG_TRACK
        stringstream _ss;
        for(auto i:globalsubgraph_.matchCandidate[1]){
         if(i==63117){
            _ss<<"first 63117"<<endl;
         }
        }
    for(auto n:globalsubgraph_.vNeighbors[63117]){
        if(n.getVertexId()==1709&&n.getMatchQueryVertexId()==0){
            _ss<<"find 1709"<<endl;
        }
    }
        Log::track1(_ss);
#endif
    return flag;
}
void CSMTopk::updateglobalVertexStarIndex(uint u1,uint v1,uint u1label,uint elabel,uint n, const std::vector<std::vector<uint>>&mcandidate) {
    globalsubgraph_.addQueryVertexCandidate(u1,v1);
    int v1size=labelToQueryVertex[query_.GetVertexLabel(u1)].size();
    globalVkMatchUk[v1].resize(n);
    for(int w=0;w<n;w++){
        globalVkMatchUk[v1][w].resize(v1size);
    }
    const std::vector<uint>&neighbors1=query_.GetNeighbors(u1);
    const std::vector<uint>&v1_nbrs=data_.GetNeighbors(v1);
    bool flagAdd;
    for(uint v1_nbr:v1_nbrs){
        for(uint n1:neighbors1){
            uint n1label=query_.GetVertexLabel(n1);
            if(std::binary_search(mcandidate[n1].begin(), mcandidate[n1].end(),v1_nbr)){
                float weight=data_.GetEdgeWeight(v1,v1_nbr);
                flagAdd=globalsubgraph_.AddEdge(u1,n1,v1,u1label,v1_nbr,n1label,elabel,weight);
                if(flagAdd) {
                    int candidate_index = queryVertexIndexInlabel[n1];
                    numupdatestar++;
                    for (int j = 0; j < query_.NumEdges(); j++) {
                        updateStarIndex( j, v1_nbr, n1, candidate_index);
                    }

                }
            }
        }
    }
    if(flagAdd) {
        int candidate_index = queryVertexIndexInlabel[u1];
        numupdatestar++;
        for (int j = 0; j < query_.NumEdges(); j++) {
            updateStarIndex( j, v1, u1, candidate_index);
        }
    }
 }
bool CSMTopk::deleteMatchRecordWithEdge(uint v1, uint v1label,uint v2, uint v2label,uint label,std::vector<int>&match) {
    bool flag = false;
    for (auto it = topKSet.begin(); it != topKSet.end();) {
        MatchRecord *record = *it;
        const std::vector<uint> &m = record->getVetex();
        bool iterflag= true;
        for (int mindex: match) {
            uint u1 = order_vs_[mindex][0];
            uint u2 = order_vs_[mindex][1];
            if ((m[u1] == v1 && m[u2] == v2) || (m[u2] == v1 && m[u1] == v2)) {
                delete record;
                record= nullptr;
                it = topKSet.erase(it);
                iterflag= false;
                flag = true;
                num_negative_results_++;
                break;
            }
        }
        if(iterflag){
            ++it;
        }
    }
    return flag;
}
bool CSMTopk::SearchMatchesWithEdge(uint m,uint v1,uint v2,float weight,uint u1,uint u2,searchType type){
    this->matchVertex(true, 0, v1, float(0));
    this->matchVertex(true, 1, v2, weight);
    if(isInsert)
        IsearchSpace+=2;
    else
        DsearchSpace+=2;
    bool isNull;
    const std::vector<uint>&uk_neighbor1=rightNeighbor[m][u1];
    total_updaterightNeighborCandidate_time.StartTimer();
    isNull=updaterightNeighborCandidate(m, u1, u2,true,v1, uk_neighbor1);
    total_updaterightNeighborCandidate_time.StopTimer();
    if(isNull)
    {
        this->popVertex(1, v2);
        this->popVertex(0, v1);
        return true;
    }
    const std::vector<uint>&uk_neighbor2=rightNeighbor[m][u2];
    total_updaterightNeighborCandidate_time.StartTimer();
    isNull=updaterightNeighborCandidate(m, u2, u1, true,v2, uk_neighbor2);
    total_updaterightNeighborCandidate_time.StopTimer();
    if(isNull)
    {
        for(int u_id:uk_neighbor1){
            int query_order_index=order_vertex_index[m][u_id];
            matchCandidate[query_order_index].clear();
        }
        this->popVertex(1, v2);
        this->popVertex(0, v1);
        return true;
    }
    searchMatches(2, m, type);
    this->popVertex( v2,m,1,uk_neighbor1);
    this->popVertex(v1,m,0,uk_neighbor2);
    return false;
 }

void CSMTopk::deleteGlobalSubgraph(uint v1, uint v2,uint elabel,float weight, std::vector<int> &match) {
    uint n=query_.NumEdges();
    bool isMatch= false;
    auto &mcandidate=globalsubgraph_.matchCandidate;
    for(auto it=match.begin();it!=match.end();it++){
        uint u1=order_vs_[*it][0];
        uint u2=order_vs_[*it][1];
        uint u1label=query_.GetVertexLabel(u1);
        uint u2label=query_.GetVertexLabel(u2);
        uint v1label=data_.GetVertexLabel(v1);
        uint v2label=data_.GetVertexLabel(v2);
        uint elabel=std::get<2>(query_.GetEdgeLabel(u1,u2));
        if(v1label!=v2label) {
            if(v1label!=u1label){
                std::swap(v1,v2);
            }
//            std::chrono::high_resolution_clock::time_point start;
            globalsubgraph_.RemoveEdge(v1,u1label,v2,u2label,u1,u2,elabel,weight);
            isMatch=deleteGlobalSubgraphHelp(*it,u1,u2,u1label,u2label,v1,v2,v1label,v2label,elabel,weight,mcandidate);
        }
        else{
            for(int i=0;i<2;i++){
                globalsubgraph_.RemoveEdge(v1,v1label,v2,v2label,u1,u2,elabel,weight);
                isMatch=deleteGlobalSubgraphHelp(*it,u1,u2,u1label,u2label,v1,v2,v1label,v2label,elabel,weight,mcandidate);
                std::swap(v1,v2);
            }
        }
       /* if(isMatch== false){
            it=match.erase(it);
        }
        else{
            ++it;
        }*/
    }
    //total_delete_update_time+= Duration2(start);
 }
 void CSMTopk::deleteUpdateglobalVertexStarIndex(uint u1,uint v1,uint v2,uint n) {
     int candidate_index= queryVertexIndexInlabel[u1];
     bool isContain;
     const std::vector<uint>&candidate=globalsubgraph_.matchCandidate[u1];
     //update u1 index
     for(int j=0;j<n;j++){
        // isContain= false;
         int vertex_index=order_vertex_index[j][u1];
         if(vertex_index==0)
             continue;
         StarGraph *s=globalStarIndex[j][vertex_index];
         if(s->getMatchDataVertexId()==v1){
             s->setStarMaxWeight(s->GetForwardNeighborNum()*mw);
             updateStarIndex(j,v1,u1,candidate_index);
             for(uint v:candidate){
                 if(globalVkMatchUk[v][j][candidate_index]>s->getStarMaxWeight()||s->getStarMaxWeight()==s->GetForwardNeighborNum()*mw)
                 {
                     s->setStarMaxWeight(globalVkMatchUk[v][j][candidate_index]);
                     s->setMatchDataVertexId(v);
                 }
             }
         }
     }
 }

 bool CSMTopk::deleteGlobalSubgraphHelp(int m,uint u1,uint u2,uint u1label,uint u2label, uint v1,uint v2,uint v1label,uint v2label,
                                          uint elabel,float weight, std::vector<std::vector<uint>>&mcandidate) {
     bool isMatch= false;
     uint n=query_.NumEdges();
     bool isContain1= true;
     bool isContain2= true;
     if(v1label!=u1label)
     {
         swap(v1,v2);
     }
     if(!std::binary_search(mcandidate[u1].begin(),mcandidate[u1].end(), v1)){
             isContain1= false;
     }
     if(!std::binary_search(mcandidate[u2].begin(),mcandidate[u2].end(), v2)){
             isContain2= false;
     }
     if(isContain1&&isContain2){
         isMatch= true;
         bool flag1= LabelFilter(v1,u1);
         bool flag2= LabelFilter(v2,u2);
         uint v1label=data_.GetVertexLabel(v1);
         uint v2label=data_.GetVertexLabel(v2);
         //globalsubgraph_.RemoveEdge(v1,v1label,v2,v2label,u1,u2,elabel,weight);
         if(flag1&&flag2){
             // still in candidate ,delete edge update bound
             deleteUpdateglobalVertexStarIndex(u1,v1,v2,n);
             deleteUpdateglobalVertexStarIndex(u2,v2,v1,n);
         }
         else if(!flag1&&flag2){
             // remove v1 update v1_nbr's index
             deleteGlobalGraphCandidateEdges(m,u1,v1,mcandidate);
         }
         else if(flag1&&!flag2){
             //remove v2 update v2_nbr's index
             deleteGlobalGraphCandidateEdges(m,u2,v2,mcandidate);
         }
         else{
             deleteGlobalGraphCandidateEdges(m,u1,v1,mcandidate);
             deleteGlobalGraphCandidateEdges(m,u2,v2,mcandidate);
         }
     }
     return isMatch;
 }



 void CSMTopk::deleteGlobalGraphCandidateEdges(uint m,uint u1,uint v1,std::vector<std::vector<uint>>&mcandidate) {
     const std::vector<uint> &neighbors1 = query_.GetNeighbors(u1);
     globalsubgraph_.deleteQueryVertexCandidate(u1, v1);
     const std::vector<uint> &v1_nbrs = data_.GetNeighbors(v1);
     bool flagDel;
     for (uint v1_nbr: v1_nbrs) {
         for (uint n1: neighbors1) {
             uint n1label = query_.GetVertexLabel(n1);
             if (std::binary_search(mcandidate[n1].begin(), mcandidate[n1].end(), v1_nbr)) {
                 uint v1label = data_.GetVertexLabel(v1);
                 uint v1_nbr_label = data_.GetVertexLabel(v1_nbr);
                 uint elabel = std::get<2>(data_.GetEdgeLabel(v1, v1_nbr));
                 uint weight = data_.GetEdgeWeight(v1, v1_nbr);
                 flagDel = globalsubgraph_.RemoveEdge(v1, v1label, v1_nbr, v1_nbr_label, u1, n1, elabel, weight);
                 if (flagDel) {
                     int candidate_index = queryVertexIndexInlabel[n1];
                     // bool isContain= false;
                     const std::vector<uint> &candidate = globalsubgraph_.matchCandidate[n1];
                     for (int j = 0; j < query_.NumEdges(); j++) {
                         //isContain= false;
                         int vertex_index = order_vertex_index[j][n1];
                         if (vertex_index == 0)
                             continue;
                         StarGraph *s = globalStarIndex[j][vertex_index];
                         if (s->getMatchDataVertexId() == n1) {
                             int vertex_index = order_vs_[j][n1];
                             if (vertex_index == 0)
                                 continue;
                             StarGraph *s = globalStarIndex[j][vertex_index];
                             s->setStarMaxWeight(s->GetForwardNeighborNum() * mw);
                             updateStarIndex(j, v1_nbr, n1, candidate_index);
                             for (uint v: candidate) {
                                 if (globalVkMatchUk[v][j][candidate_index] > s->getStarMaxWeight() ||
                                     s->getStarMaxWeight() == s->GetForwardNeighborNum() * mw) {
                                     s->setStarMaxWeight(globalVkMatchUk[v][j][candidate_index]);
                                 }
                             }

                         }

                     }
                 }
             }
         }
     }
 }
