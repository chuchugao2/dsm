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
#include "../utils/Timer.h"

class Graphflow : public matching
{
public:
    // a list of matching orders starting from each query edge
    // the first matching order also applies to the initial matching
    std::vector<std::vector<uint> > order_vs_; //匹配的点顺序
    std::vector<std::vector<uint> > order_csrs_;//匹配节点前向邻居
    std::vector<std::vector<uint> > order_offs_;//匹配节点的索引
    std::vector<std::vector<uint>> order_vertex_index;//每个节点在该匹配序中的位置 order_vertex_index[0][u1]表示第一个匹配序中，u1的索引位置
//        std::unordered_map<std::pair<uint,uint>,std::vector<MatchRecord*>,pair_hash> edgeMaps;//对应于每条边的辅助索引结构
    // uint s,t;//findMatch时更新键值
    //记录top k记录集合
    std::vector<MatchRecord*> topKSet;
    //记录所有匹配的结果
    std::vector<std::vector<StarGraph*>>globalStarIndex;//记录每种匹配下每个节点左邻居以及最大权值
    std::vector<SingleCandidate>match;//每个节点的匹配结果  vertex density
    std::vector<std::vector<SingleCandidate>>matchCandidate;//id density
    std::vector<float>suffixMax;
    std::vector<float>isolatedMax;
    std::vector<std::vector<std::vector<uint>>>rightNeighbor;//匹配索引号，id号
    //std::vector<LocalIndex>queryLocalIndexs;
    std::vector<std::vector<std::vector<int>>>globalVkMatchUk;//<vk,ak,uk>
    std::vector<std::vector<uint>>labelToQueryVertex;//每个标签对应的查询点标签
    std::vector<uint>queryVertexIndexInlabel;//每个查询点在label数组中的索引号
    std::vector<float>LocalStarIndex;//局部索引
    std::vector<std::vector<int>>matchLeftNeighborSum;//所有节点左邻居的个数
    std::vector<std::vector<size_t>>matchVetexLeftNeighbor;//所有匹配序列中左邻居组合数
    std::vector<vector<StarGraph*>>matchVetexSumweight;//每种组合更新得到的最大权值
    std::vector<std::vector<size_t>>leftNeighborIdSum;//每个节点左邻居id和
    bool isUpdateIntopkset;
    int numAddTopk=0;
    int allMatchFind=0;
    int sumAllMatchFind=0;
    int sumDeleteallMatchFind=0;
    int numupdatestar=0;
    long IsearchSpace=0,DsearchSpace=0,IdeterminCandite=0,DdeterminCandite=0;

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
    void InitialTopK(const std::string &path) override;//得到初始化之后的Top k结果集合
    void updateTopK() override;
    void deleteEdge(uint v1,uint v2) override;
    void deleteUpdateTopK() override;
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
    void PrintAverageTime(int len);
private:
    void GenerateMatchingOrder();
    void FindMatches(uint flag,uint order_index, uint depth,
                     std::vector<uint> m, size_t &num_results,float density_s); //flag==0 initial flag==1 update
    int addMatchRecords(MatchRecord* r);//1 表示成功插入 2表示节点重复 3表示节点已比第kth小
    void addStarGraph(StarGraph *s);
    void CreateStarIndex();
    float GetBackWeight(uint order_index,uint depth);
    void updateStarIndex(uint match_index,uint caddidate_v,const std::vector<uint>&canditeQueryVertexs);
    void updateStarIndex(uint match_index,uint caddidate_v,uint candidate_u,int candidate_v_index);
    vector<int> EdgeisInMatchOrder(Edge *edge);
    vector<int> EdgeisInMatchOrder(uint v1,uint v2,uint v1label,uint v2label,uint velabel);
    void searchMatches(int depth,uint matchorderindex,searchType flag);
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
   void createGlobalSubgraph();//构建全局子图
   bool updateGlobalSubgraph(uint v1,uint v2,uint label,float weight, std::vector<int>&match);
   bool updateGlobalGraphHelp(int m,uint u1,uint u2,uint u1label,uint u2label, uint v1,uint v2,uint v1label,uint v2label,uint elabel,const std::vector<std::vector<uint>>&mcandidate,bool &flag);
   void updateglobalVertexStarIndex(uint u1,uint v1,uint u1label,uint elabel,uint n, const std::vector<std::vector<uint>>&mcandidate);//新增v1候选节点，并更新全局索引
   void deleteGlobalSubgraph(uint v1, uint v2,uint elabel,float weight, std::vector<int> &match);
   void deleteUpdateglobalVertexStarIndex(uint u1,uint v1,uint v2,uint n);
   bool deleteGlobalSubgraphHelp(int m,uint u1,uint u2,uint u1label,uint u2label, uint v1,uint v2,uint v1label,uint v2label,
                                 uint elabel,float weight, std::vector<std::vector<uint>>&mcandidate);
   void deleteGlobalGraphCandidateEdges(uint m,uint u1,uint v1,std::vector<std::vector<uint>>&mcandidate);
   bool deleteMatchRecordWithEdge(uint v1, uint v1label,uint v2, uint v2label,uint label,std::vector<int> &match);
   bool SearchMatchesWithEdge(uint m,uint v1,uint v2,float weight,uint u1,uint u2,searchType type);

};

#endif //MATCHING_GRAPHFLOW
