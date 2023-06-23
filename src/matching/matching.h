#ifndef CSM_MATCHING_H
#define CSM_MATCHING_H

#include <vector>

#include "../utils/types.h"
#include "../graph/graph.h"
#include "../graph/Subgraph.h"


class matching
{
protected:
    Graph& query_; //��ѯͼ
    Graph& data_;//����ͼ
    Subgraph& globalsubgraph_;//ȫ�ֺ�ѡ����


    // config
    const size_t max_num_results_; //��ÿ�θ��µ����Ľ������
    const bool print_preprocessing_results_; //�Ƿ��ӡԤ������
    const bool print_enumeration_results_;//�Ƿ��ӡƥ����
    const bool homomorphism_;//�Ƿ�����̬ͬ

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

    virtual void Preprocessing();//Ԥ����
    virtual void InitialMatching(const std::string &path);//��ʼ��ƥ��˳��
    virtual void AddEdge(uint v1, uint v2, uint label,float weight);//���ӱ�
    virtual void RemoveEdge(uint v1, uint v2,uint label);//ȥ����
    virtual void AddVertex(uint id, uint label);//���ӽڵ�
    virtual void RemoveVertex(uint id);//ȥ���ڵ�

    virtual void GetMemoryCost(size_t &num_edges, size_t &num_vertices);//�õ��ڴ滨��


    virtual void clearAllMatches();
    virtual void InitialTopK(const std::string &path);
    virtual void updateTopK(uint num);//����Top k����
    virtual void deleteEdge(uint v1,uint v2);
    virtual void deleteUpdateTopK();
    virtual void PrintAverageTime(int len);
    // get execution info
    void GetNumInitialResults(size_t &num_initial_results);//�õ���ʼ��������
    void GetNumPositiveResults(size_t &num_positive_results);//��ƥ�������
    void GetNumNegativeResults(size_t &num_negative_results);//��ƥ�������
    void clearPositiveNum(); //��� ��ƥ�����
    void PrintCounter();//��ӡʣ���ִ����Ϣ


};

#endif //CSM_MATCHING_H
