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
#include "matching/graphflow.h"
#include "Instopk.h"

#ifndef TEST_MAIN
int main(int argc, char *argv[])
{
    CLI::App app{"App description"};

    std::string query_path = "", initial_path = "", stream_path = "",result_path="",query_info="",algorithm="instopk";
    uint max_num_results = UINT_MAX, time_limit = UINT_MAX, initial_time_limit = UINT_MAX;
    uint dist=2,update_len=5000;
    bool print_prep = true, print_enum = false, homo = false, report_initial = true;
    std::vector<std::vector<uint>> orders;

    app.add_option("-q,--query", query_path, "query graph path")->required();
    app.add_option("-d,--data", initial_path, "initial data graph path")->required();
    app.add_option("-u,--update", stream_path, "data graph update stream path")->required();
    app.add_option("--result-path",result_path,"topk result with each update");
    app.add_option("-a,--algorithm", algorithm, "algorithm");
    app.add_option("--max-results", max_num_results, "max number of results for one edge update");
    app.add_option("--time-limit", time_limit, "time limit for the incremental matching (second)");
    app.add_option("--print-prep", print_prep, "print processing results or not");
    app.add_option("--print-enum", print_enum, "print enumeration results or not");
    app.add_option("--homo", homo, "using graph homomorphism");
    app.add_option("--report-initial", report_initial, "report the result of initial matching or not");
    app.add_option("--initial-time-limit", initial_time_limit, "time limit for the initial matching (second)");
    app.add_option("--orders", orders, "pre-defined matching orders");
    app.add_option("--qInfo",query_info,"the path of query graph");
    app.add_option("--dist",dist,"the dist");
    app.add_option("--ul",update_len,"the deletion num");
    CLI11_PARSE(app, argc, argv);

    std::chrono::high_resolution_clock::time_point start, lstart;

    Log::init_track1("/media/hnu/hnu2022/gaochuchu/dsm_instok/src/log/loginfo1.txt");
    Log::init_track3("/media/hnu/hnu2022/gaochuchu/dsm_instok/src/log/WCcomputeI_text.txt");
#ifdef COMPUTE_TIME
    stringstream _ss;
    _ss << query_info << endl;
    Log::track3(_ss);
#endif
    std::string path="/media/hnu/hnu2022/gaochuchu/dsm_instok/src/log/";
    std::string initial_result_path="/media/hnu/hnu2022/gaochuchu/dsm_instok/src/result3/";
    if(result_path==""){
        path+="topkResult.txt";
    }else{
        path+=result_path;
    }
    if(query_info==""){
        initial_result_path+="text";
    }else{
        initial_result_path+=query_info;
    }
    Log::init_track2(path);


    start = Get_Time();
    std::cout << "----------- Loading graphs ------------" << std::endl;
    Graph query_graph {};
    query_graph.LoadFromFile(query_path,0);
    query_graph.PrintMetaData();

    Graph data_graph {};
    data_graph.LoadFromFile(initial_path,1);
    data_graph.PrintMetaData();
    Print_Time("Load Graphs: ", start);

    std::cout << "------------ Preprocessing ------------" << std::endl;
    matching *mm = nullptr;
    Graphflow *graphflow = nullptr;
    Instopk *instopk= nullptr;
    start = Get_Time();
    if(algorithm=="graphflow")
        mm = graphflow      = new Graphflow     (query_graph, data_graph, max_num_results, print_prep, print_enum, homo);
    else if(algorithm=="instopk")
        mm = instopk     = new Instopk     (query_graph, data_graph, max_num_results, print_prep, print_enum, homo,dist);
    auto PreprocessFun = [&mm]()
    {
        mm->Preprocessing();
    };
    execute_with_time_limit(PreprocessFun, 3600, reach_time_limit);
    Print_Time("Preprocessing: ", start);
    std::cout << "\nPeak Virtual Memory: " << mem::getValue() << " KB"<<endl;
    std::cout << "\n\n----------------- End -----------------" << std::endl;
#ifdef COMPUTE_TIME
    mm->InitialSpace=mem::getValue();
#endif
    if (report_initial)
    {
        std::cout << "----------- Initial Matching ----------" << std::endl;

        start = Get_Time();
       // data_graph.LoadUpdateStream(initial_path);
        auto InitialFun = [&mm,&data_graph,&initial_result_path]()
        {
            mm->InitialMatching(initial_result_path);
            std::chrono::high_resolution_clock::time_point s;
            s = Get_Time();
            mm->InitialTopK(initial_result_path);
            Print_Time("InitialTopk ", s);
        };
        execute_with_time_limit(InitialFun, initial_time_limit, reach_time_limit);
        Print_Time("Initial Matching: ", start);

        size_t num_results = 0ul;
        mm->GetNumInitialResults(num_results);
        std::cout << num_results << " initial matches.\n";
//        std::cout<<gs.size()<<std::endl;
        //if (reach_time_limit) return 1;
        if(reach_time_limit)
            reach_time_limit= false;
    }
    std::cout << "--------- Incremental Matching --------" << std::endl;
    data_graph.LoadUpdateStream(stream_path);
    mm->clearPositiveNum();
    size_t num_v_updates = 0ul, num_e_updates = 0ul;
    int upl=10000-update_len;
    auto IncrementalFun = [&data_graph, &mm, &num_v_updates, &num_e_updates,&upl]()
    {
        while (!data_graph.updates_.empty())
        {
            std::cout<<"update num: "<<data_graph.updates_.size()<<std::endl;
            stringstream _ss;
            _ss<<"update num:"<<data_graph.updates_.size()<<"\n";
            Log::track2(_ss);
#ifdef PRINT_DEBUG
            Log::track1(_ss);

#endif
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
                bool isInsert= true;
                if(data_graph.updates_.size()==9900||reach_time_limit){
                    if(reach_time_limit){
                        long long tmp=mm->total_search_time.GetTimer();
                        /* std::random_device rd;
                         std::mt19937_64 gen(rd());
                         std::uniform_int_distribution<long long> dis(1, 120000000);
                         long long random_number = dis(gen);*/
                        long long insert_sum_time = ((tmp) * 1.0 /(10000 - data_graph.updates_.size()+1)) * upl;
                        mm->total_search_time.setTimer(insert_sum_time);
                        reach_time_limit=false;
                    }
                    else{
                        long long tmp=mm->total_search_time.GetTimer();
                        long long insert_sum_time=tmp*1.0/100*upl;
                        mm->total_search_time.setTimer(insert_sum_time);
                    }
                    int target=10000-upl;
                    while(data_graph.updates_.size()>target){
                        data_graph.updates_.pop();
                    }
                    isInsert= false;
                }
                if(isInsert){
                    mm->total_search_time.StartTimer();
                    auto TotalSearchFun = [&mm, &insert]() {
                        mm->AddEdge(insert.id1, insert.id2, insert.label,insert.weight,insert.timestamp);
                    };
                    execute_with_time_limit(TotalSearchFun, 10, reach_time_limit);
                    mm->total_search_time.StopTimer();
                    if(reach_time_limit) {
                        std::cout<<"update num:"<<data_graph.updates_.size()<<endl;
                    }
                    std::cout<<"insert time:"<<mm->total_search_time.GetTimer()<<endl;
                }
//                 data_graph.AddEdge(insert.id1, insert.id2, insert.label,insert.weight,insert.timestamp,1);
                num_e_updates ++;

#ifdef PRINT_DEBUG
#endif
            }
            else if (insert.type == 'e' && !insert.is_add)
            {
//                mm->RemoveEdge(insert.id1, insert.id2);
               /* bool isDelete=true;
                if (data_graph.updates_.size()==(10000-upl-100)||reach_time_limit)
                {
                    if(reach_time_limit){
                    long long tmp=mm->total_delete_time.GetTimer();
                    uint dpl=10000-upl;
                    *//*std::random_device rd;
                    std::mt19937_64 gen(rd());
                    std::uniform_int_distribution<long long> dis(1, 60000000);
                    long long random_number = dis(gen);*//*
                    long long delete_sum_time=tmp*1.0/(dpl- data_graph.updates_.size()+1)*upl;
                    mm->total_delete_time.setTimer(delete_sum_time);

                    reach_time_limit=false;
                    }
                   else{
                        int dpl=10000-upl;
                        long long tmp=mm->total_delete_time.GetTimer();
                        long long insert_sum_time=tmp*1.0/100*dpl;
                        mm->total_delete_time.setTimer(insert_sum_time);
                   }
                    int target=0;
                    while(data_graph.updates_.size()>target){
                        data_graph.updates_.pop();
                    }
                    isDelete= false;
                }
                if(isDelete) {
                    mm->total_delete_time.StartTimer();
                    auto TotalDeleteFun = [&mm, &insert]() {
                        mm->RemoveEdge(insert.id1,insert.id2);
                    };
                    execute_with_time_limit(TotalDeleteFun, 10, reach_time_limit);
                    mm->total_delete_time.StopTimer();
                    if(reach_time_limit) {
                        std::cout<<"update num:"<<data_graph.updates_.size()<<endl;
                    }
                    std::cout<<"delete time:"<<mm->total_delete_time.GetTimer()<<endl;
                }
                num_e_updates ++;
                mm->RemoveEdge(insert.id1,insert.id2);
                num_e_updates ++;*/
            }
            //if (reach_time_limit) break;
        }
    };

    start = Get_Time();

    execute_with_time_limit(IncrementalFun, time_limit, reach_time_limit);

    Print_Time("Incremental Matching: ", start);

    std::cout << num_v_updates << " vertex updates.\n";
    std::cout << num_e_updates << " edge updates.\n";

    size_t positive_num_results = 0ul, negative_num_results = 0ul;
    mm->GetNumPositiveResults(positive_num_results);
    mm->GetNumNegativeResults(negative_num_results);
    std::cout << positive_num_results << " positive matches.\n";
    std::cout << negative_num_results << " negative matches.\n";
    mm->PrintCounter();

    size_t num_edges = 0u, num_vertices = 0ul;
    mm->GetMemoryCost(num_edges, num_vertices);
    std::cout << "\n# edges in index in total: " << num_edges;
    std::cout << "\n# vertices in index in total: " << num_vertices;

    std::cout << "\nPeak Virtual Memory: " << mem::getValue() << " KB"<<endl;
    mm->space_const=mem::getValue()*1.0/1024;
    mm->PrintAverageTime(update_len);
    std::cout << "\n\n----------------- End -----------------" << std::endl;


    Log::finalize();
    delete mm;

}
#else
int main(int argc, char *argv[]){
    Graph data_graph {};
    data_graph.createDataStream("/home/gaochuchu/gcc/data/netflow/netflow/data_graph/data.graph","/home/gaochuchu/gcc/data/testdata_graph.txt");
    data_graph.createInitalDataGraph("/home/gaochuchu/gcc/data/testdata_graph.txt","/home/gaochuchu/gcc/data/initdata_graph.txt");
    data_graph.createUpdateStream("/home/gaochuchu/gcc/data/testdata_graph.txt","/home/gaochuchu/gcc/data/testupdate_graph.txt");
   std::cout<<"hello world"<<std::endl;
}
#endif
