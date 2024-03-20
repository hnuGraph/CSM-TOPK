#include <chrono>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include "utils/CLI11.hpp"
#include "utils/globals.h"
#include "utils/types.h"
#include "utils/Log.h"
#include "graph/graph.h"
#include "matching/matching.h"
#include "matching/CSMTopk.h"

int main(int argc, char *argv[])
{
    CLI::App app{"App description"};

    std::string query_path = "", initial_path = "", stream_path = "",query_info="";
    uint max_num_results = UINT_MAX, time_limit = UINT_MAX, initial_time_limit = UINT_MAX;
    uint update_len=5000,k_size=100;
    bool print_prep = true, print_enum = false, homo = false, report_initial = true;
    if(print_enum)
        print_result= true;

    std::string initial_result_path="";
//    Log::init_track1("");
//    Log::init_track3("");
    app.add_option("-q,--query", query_path, "query graph path")->required();
    app.add_option("-d,--data", initial_path, "initial data graph path")->required();
    app.add_option("-u,--update", stream_path, "data graph update stream path")->required();
    app.add_option("--max-results", max_num_results, "max number of results for one edge update");
    app.add_option("--time-limit", time_limit, "time limit for the incremental matching (second)");
    app.add_option("--print-prep", print_prep, "print processing results or not");
    app.add_option("--print-enum", print_enum, "print enumeration results or not");
    app.add_option("--print-result", print_result, "print update top k  results or not");
    app.add_option("--report-initial", report_initial, "report the result of initial matching or not");
    app.add_option("--initial-time-limit", initial_time_limit, "time limit for the initial matching (second)");
    app.add_option("--qInfo",query_info,"the path of query graph");
    app.add_option("--ul",update_len,"the deletion len");
    app.add_option("--ksize",k_size,"the k size of query");

    CLI11_PARSE(app, argc, argv);
    k=k_size;
#ifdef COMPUTE_TRACK
    stringstream _ss;
    _ss<<query_info<<endl;
    Log::track3(_ss);
#endif
    std::chrono::high_resolution_clock::time_point start;
    start = Get_Time();
    std::cout << "----------- Loading graphs ------------" << std::endl;
    Graph query_graph {};
    query_graph.LoadFromFile(query_path,0);
    query_graph.PrintMetaData();

    Graph data_graph {};
    data_graph.LoadFromFile(initial_path,1);
    data_graph.PrintMetaData();
    Print_Time("Load Graphs: ", start);
    Subgraph globalsubgraph(data_graph.NumVertices(),query_graph.NumVertices());
    std::cout << "------------ Preprocessing ------------" << std::endl;
    matching *mm = nullptr;
    CSMTopk *csmtopk = nullptr;
    start = Get_Time();
        mm = csmtopk      = new CSMTopk    (query_graph, data_graph,globalsubgraph, max_num_results, print_prep, print_enum, homo);
    mm->Preprocessing();
    Print_Time("Preprocessing1: ", start);
//    stringstream _ss1;
//    _ss1<<mm->total_index_time.GetTimer()<<endl;
//    Log::track3(_ss1);

    if (report_initial)
    {
        std::cout << "----------- Initial Matching ----------" << std::endl;

        start = Get_Time();
       // data_graph.LoadUpdateStream(initial_path);
        auto InitialFun = [&mm,&data_graph,&initial_result_path]()
        {
            mm->InitialMatching(initial_result_path);
            mm->InitialTopK(initial_result_path);
        };
        execute_with_time_limit(InitialFun, initial_time_limit, reach_time_limit);
        Print_Time("Initial Matching: ", start);

        size_t num_results = 0ul;
        mm->GetNumInitialResults(num_results);
        std::cout << num_results << " initial matches.\n";
        if (reach_time_limit) return 1;
    }
    std::cout << "--------- Incremental Matching --------" << std::endl;
    data_graph.LoadUpdateStream(stream_path);
    uint upl=10000-update_len;
    mm->clearPositiveNum();
    size_t num_v_updates = 0ul, num_e_updates = 0ul;

    auto IncrementalFun = [&data_graph, &mm, &num_v_updates, &num_e_updates,&upl]()
    {
        while (!data_graph.updates_.empty())
        {
                if(data_graph.updates_.size()==upl){
                    mm->isInsert= false;
                    mm->Itotal_updaterightNeighborCandidate_time=mm->total_updaterightNeighborCandidate_time.GetTimer();
                    mm->Itotal_densityfilter_time=mm->total_densityFilter_time.GetTimer();
                    mm->total_updaterightNeighborCandidate_time.clearTimer();
                    mm->total_densityFilter_time.clearTimer();
                }
if(print_result)
            std::cout<<"update num: "<<data_graph.updates_.size()<<std::endl;

            InsertUnit insert = data_graph.updates_.front();
            data_graph.updates_.pop();

            if (insert.type == 'v' && insert.is_add)
            {
                mm->AddVertex(insert.id1, insert.label);
                num_v_updates ++;
            }
            else if (insert.type == 'v' && !insert.is_add)
            {
                mm->RemoveVertex(insert.id1);
                num_v_updates ++;

            }
            else if (insert.type == 'e' && insert.is_add)
            {
                mm->AddEdge(insert.id1, insert.id2, insert.label,insert.weight);
                num_e_updates ++;

            }
            else if (insert.type == 'e' && !insert.is_add)
            {
                mm->total_delete_all.StartTimer();
                mm->RemoveEdge(insert.id1,insert.id2,insert.label);
                mm->total_delete_all.StopTimer();
                num_e_updates ++;
            }
            if (reach_time_limit) break;
        }
    };

    start = Get_Time();
    execute_with_time_limit(IncrementalFun, time_limit, reach_time_limit);
    Print_Time("Incremental Matching: ", start);
    std::cout << num_e_updates << " edge updates.\n";
    size_t positive_num_results = 0ul, negative_num_results = 0ul;
    size_t num_edges = 0u, num_vertices = 0ul;
    mm->space_cost=mem::getValue();
    mm->PrintAverageTime(update_len);
    std::cout << "\n\n----------------- End -----------------" << std::endl;

    Log::finalize();
    delete mm;

}
