#ifndef CSM_MATCHING_H
#define CSM_MATCHING_H

#include <vector>

#include "../utils/types.h"
#include "../graph/graph.h"
#include "../utils/Timer.h"


class matching {
protected:
    Graph &query_; //查询图
    Graph &data_;//数据图


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
    matching(Graph &query_graph, Graph &data_graph,
             size_t max_num_results = ULONG_MAX,
             bool print_preprocessing_results = true,
             bool print_enumeration_results = false,
             bool homomorphism = false);

    virtual ~matching() = default;

    virtual void Preprocessing();//预处理
    virtual void InitialMatching(const std::string &path);//初始化匹配顺序

    virtual void AddEdge(uint v1, uint v2, uint label, float weight, uint timestamp);//增加边
    virtual void AddEdgeWithGlobalIndex(uint v1, uint v2, uint label, float weight, uint timestamp);

    virtual void RemoveEdge(uint v1, uint v2, uint label);//去除边
    virtual void AddVertex(uint id, uint label);//增加节点
    virtual void RemoveVertex(uint id);//去除节点

    virtual void GetMemoryCost(size_t &num_edges, size_t &num_vertices);//得到内存花费


    virtual void clearAllMatches();

    virtual void InitialTopK(const std::string &path);

    virtual void updateTopK();//更新Top k序列
    virtual void deleteEdge(uint v1, uint v2);

    virtual void deleteUpdateTopK();

    virtual void PrintAverageTime(int len);

    // get execution info
    void GetNumInitialResults(size_t &num_initial_results);//得到初始化的数量
    void GetNumPositiveResults(size_t &num_positive_results);//正匹配的数量
    void GetNumNegativeResults(size_t &num_negative_results);//负匹配的数量
    void clearPositiveNum(); //清空 正匹配个数
    void PrintCounter();//打印剩余的执行信息
    Timer total_search_time, total_print_time, total_densityFilter_time, total_update_globalIndex_time, total_updaterightNeighborCandidate_time,
            total_delete_time, total_delete_update_time,total_test1,total_test2,total_test3,total_test4,total_test5,total_test6, total_test7;
    long long Itotal_densityfilter_time=0, Itotal_updaterightNeighborCandidate_time=0,Itotal_test1=0,Itotal_test2=0,
    Itotal_test3=0,Itotal_test4=0,Itotal_test5=0,Itotal_test6=0,Itotal_test7=0;
    bool isInsert= true;


};

#endif //CSM_MATCHING_H
