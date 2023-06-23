#ifndef MATCHING_GRAPHFLOW
#define MATCHING_GRAPHFLOW

#include <vector>
#include "../utils/types.h"
#include "../graph/graph.h"
#include "../matching/matching.h"
#include "unordered_map"
#include "sstream"
#include "../utils/Log.h"
#include "../utils/globals.h"
#include "../graph/MatchRecord.h"
#include "../graph/StarGraph.h"
#include "algorithm"
#include "cfloat"
#include "LocalIndex.h"
#include "SingleCandidate.h"


class Graphflow : public matching
{
public:
    // a list of matching orders starting from each query edge
    // the first matching order also applies to the initial matching
    std::vector<std::vector<uint> > order_vs_; //ƥ��ĵ�˳��
    std::vector<std::vector<uint> > order_csrs_;//ƥ��ڵ�ǰ���ھ�
    std::vector<std::vector<uint> > order_offs_;//ƥ��ڵ������
    std::vector<std::vector<uint>> order_vertex_index;//ÿ���ڵ��ڸ�ƥ�����е�λ�� order_vertex_index[0][u1]��ʾ��һ��ƥ�����У�u1������λ��
//        std::unordered_map<std::pair<uint,uint>,std::vector<MatchRecord*>,pair_hash> edgeMaps;//��Ӧ��ÿ���ߵĸ��������ṹ
    // uint s,t;//findMatchʱ���¼�ֵ
    //��¼top k��¼����
    std::vector<MatchRecord*> topKSet;
    //��¼����ƥ��Ľ��
    std::vector<std::vector<StarGraph*>>globalStarIndex;//��¼ÿ��ƥ����ÿ���ڵ����ھ��Լ����Ȩֵ
    std::vector<SingleCandidate>match;//ÿ���ڵ��ƥ����  vertex density
    std::vector<std::vector<SingleCandidate>>matchCandidate;//id density
    std::vector<float>suffixMax;
    std::vector<float>isolatedMax;
    std::vector<std::vector<std::vector<uint>>>rightNeighbor;//ƥ�������ţ�id��
    //std::vector<LocalIndex>queryLocalIndexs;
    std::vector<std::vector<std::vector<int>>>globalVkMatchUk;//<vk,ak,uk>
    std::vector<std::vector<uint>>labelToQueryVertex;//ÿ����ǩ��Ӧ�Ĳ�ѯ���ǩ
    std::vector<uint>queryVertexIndexInlabel;//ÿ����ѯ����label�����е�������
    std::vector<float>LocalStarIndex;//�ֲ�����
    std::vector<std::vector<int>>matchLeftNeighborSum;//���нڵ����ھӵĸ���
    std::vector<std::vector<int>>matchVetexLeftNeighbor;//����ƥ�����������ھ������
    std::vector<vector<float>>matchVetexSumweight;//ÿ����ϸ��µõ������Ȩֵ
    std::vector<std::vector<int>>leftNeighborIdSum;//ÿ���ڵ����ھ�id��
    long long total_search_time=0;
    long long total_update_gloabalsubgraph_time=0;
    long long total_print_time=0;
    long long total_densityFilter_time=0;
    long long total_update_globalIndex_time=0;
    long long total_update_localIndex_time=0;
    long long total_addMatchResult_time=0;
    long long subtime_in_local_index=0;


public:
    Graphflow(Graph& query_graph, Graph& data_grasph,Subgraph& global_subgraph, uint max_num_results,
              bool print_prep, bool print_enum, bool homo);
    ~Graphflow() override ;

    void Preprocessing() override;
    void InitialMatching(const std::string &path) override;

    void AddEdge(uint v1, uint v2, uint label,float weight) override;
    void RemoveEdge(uint v1, uint v2,uint label) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    void InitialTopK(const std::string &path) override;//�õ���ʼ��֮���Top k�������
    void updateTopK(uint num) override;
    void deleteEdge(uint v1,uint v2) override;
    void deleteUpdateTopK() override;
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
    void PrintAverageTime(int len);
private:
    void GenerateMatchingOrder();
    void FindMatches(uint flag,uint order_index, uint depth,
                     std::vector<uint> m, size_t &num_results,float density_s); //flag==0 initial flag==1 update
    int addMatchRecords(MatchRecord* r);//1 ��ʾ�ɹ����� 2��ʾ�ڵ��ظ� 3��ʾ�ڵ��ѱȵ�kthС
    void addStarGraph(StarGraph *s);
    void CreateStarIndex();
    float GetBackWeight(uint order_index,uint depth);
    void updateStarIndex(uint match_index,uint caddidate_v,const std::vector<uint>&canditeQueryVertexs);
    void updateStarIndex(bool isAdd,uint match_index,uint caddidate_v,uint candidate_u,int candidate_v_index);
    vector<int> EdgeisInMatchOrder(Edge *edge);
    vector<int> EdgeisInMatchOrder(uint v1,uint v2,uint v1label,uint v2label,uint velabel);
    void searchMatches(int depth,uint matchorderindex,searchType flag,float maxWeight);
    bool LabelFilter(uint data_v,uint query_v);
    void matchVertex(bool isFirstEdge,uint depth,uint data_v,float w);
    void matchVertex(int depth);
    void popVertex(uint depth,uint data_v);
    void popVertex(uint data_v,uint matchorderindex,uint depth, const std::vector<uint>&uk_neighbor);
    void popVertex(uint depth);
    void densityFilter(uint matchorder_index,uint depth, std::vector<SingleCandidate>&singleVertexCandidate);
    void combination_helper(std::vector<std::vector<int>>& result, std::vector<int>& current, const std::vector<std::vector<uint>>& nums, int k);
    std::vector<std::vector<int>> combination(const std::vector<std::vector<uint>>& nums);
    bool LDVertexCandidateCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection,  std::vector<uint> &intersectresult);
    void setSingleVertexByIntersectionResult( std::vector<tuple<int,int,float>>&singleVertex,std::vector<uint> &intersectresult,std::vector<int>&r);
    void setIsolateVertexMatchResult(std::vector<int>&r,std::vector<int>&isolateVertex,float density);
    void setBatchVisited(std::vector<int>&r,bool flag);
    void recoverIsolateVertexMatchResult(std::vector<int>&IsolateVertexs);
    void sychronizeSingleVertexAndCandidate( std::vector<tuple<int,int,float>>&singleVertex,std::vector<uint> &intersectresult);
    void addMatchResult(uint matchorderindex,searchType type);
    std::vector<tuple<std::vector<int>,int,float>>combinationMatchResult(std::vector<std::vector<tuple<int,int,float>>>combinezIsolateVertexs);
    void combinationMatchResultHelp(std::vector<tuple<std::vector<int>,int,float>>&result,std::vector<int>&current,
                                    std::vector<std::vector<tuple<int,int,float>>>&combinezIsolateVertexs,int k,int tmin,float density
    );
    float findWeightBeforeIsolated();
    void CatesianProductWithIndex(int matchorderindex,searchType type,int curIndex,int depth,int len,int*hash,std::vector<std::vector<SingleCandidate>>&combinezIsolateVertexs,std::vector<int>&isolateVertexs,float &weight);
    int  findTboundMaxIndex(float *Tbound,int*hash,int*nocan,std::vector<std::vector<SingleCandidate>>&combinezIsolateVertexs,int len);
    bool isnoNextVertex(int*noscan,int len);
    void addMatchResultWithHeap(uint matchorderindex,searchType type);
    void CatesianProductWithHeap(int matchorderindex, searchType type, int depth, int len, int *hash,
                                 std::vector<std::vector<SingleCandidate>> &combinezIsolateVertexs,
                                 std::vector<int> &isolateVertexs, std::vector<int> &isolatedIndex,float &weight);
    void createLabelToQueryVertex();

   bool updaterightNeighborCandidate(int matchorderindex,uint uk,uint uk_neigh,bool isFirstEdge,uint vk,const std::vector<uint>&uk_neighbor);
   void InitialLocalIndex(int matchorderindex);
   void getIntersetSingleCandidate( std::vector<SingleCandidate>&candidates,int matchorderindex,int depth);
   void createGlobalSubgraph();//����ȫ����ͼ
   bool updateGlobalSubgraph(uint v1,uint v2,uint label,float weight, std::vector<int>&match);
   bool updateGlobalGraphHelp(int m,uint u1,uint u2,uint u1label,uint u2label, uint v1,uint v2,uint v1label,uint v2label,uint elabel,const std::vector<std::vector<uint>>&mcandidate,bool &flag);
   void deleteGlobalSubgrah(uint v1,uint v2,std::vector<uint>&match);
   //bool deleteGlobalGraphHelp(uint u1,uint u2,uint u1label,uint u2label, uint v1,uint v2,uint v1label,uint v2label,uint elabel,const std::vector<std::vector<uint>>&mcandidate,bool &flag);
    bool deleteMatchRecordWithEdge(uint v1, uint v1label,uint v2, uint v2label,uint label, float &maxWeight);

};

#endif //MATCHING_GRAPHFLOW
