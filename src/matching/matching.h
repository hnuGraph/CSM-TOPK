#ifndef CSM_MATCHING_H
#define CSM_MATCHING_H

#include <vector>

#include "../utils/types.h"
#include "../graph/graph.h"
#include "../graph/Subgraph.h"
#include "../utils/Timer.h"


class matching
{
protected:
    Graph& query_; //查询图
    Graph& data_;//数据图
    Subgraph& globalsubgraph_;//全局候选索引


    // config
    const size_t max_num_results_; //在每次更新的最大的结果数量
    const bool print_preprocessing_results_; //是否打印预处理结果
    const bool print_enumeration_results_;//是否打印匹配结果
    const bool homomorphism_;//是否允许同态

    // execution info
    std::vector<bool> visited_;
    size_t num_initial_results_;
    size_t num_positive_results_;
    size_t num_negative_results_;
    size_t num_intermediate_results_before_index_check_;
    size_t num_intermediate_results_after_index_check_;
    size_t num_intermediate_results_after_joinability_check_;
    size_t num_intermediate_results_after_visit_check_;
    size_t num_intermediate_results_with_empty_candidate_set_;
    size_t num_intermediate_results_without_results_;


public:
    matching(Graph& query_graph, Graph& data_graph,Subgraph& global_subgraph,
             size_t max_num_results = ULONG_MAX,
             bool print_preprocessing_results = true,
             bool print_enumeration_results = false,
             bool homomorphism = false);
    virtual ~matching() = default;

    virtual void Preprocessing();
    virtual void InitialMatching(const std::string &path);//Initial matching
    virtual void AddEdge(uint v1, uint v2, uint label,float weight);
    virtual void RemoveEdge(uint v1, uint v2,uint label);
    virtual void AddVertex(uint id, uint label);
    virtual void RemoveVertex(uint id);
    virtual void GetMemoryCost(size_t &num_edges, size_t &num_vertices);
    virtual void clearAllMatches();
    virtual void InitialTopK(const std::string &path);
    virtual void updateTopK();
    virtual void deleteUpdateTopK();
    virtual void PrintAverageTime(int len);
    // get execution info
    void GetNumInitialResults(size_t &num_initial_results);
    void GetNumPositiveResults(size_t &num_positive_results);
    void GetNumNegativeResults(size_t &num_negative_results);
    void clearPositiveNum();
    void PrintCounter();
    Timer total_search_time, total_print_time, total_densityFilter_time, total_update_globalIndex_time, total_updaterightNeighborCandidate_time,
            total_delete_time, total_delete_update_time,total_delete_all,total_index_time;
    long long Itotal_densityfilter_time=0, Itotal_updaterightNeighborCandidate_time=0;
    bool isInsert= true;
    float space_cost=0;
    int curupdatenum=0;

};

#endif //CSM_MATCHING_H
