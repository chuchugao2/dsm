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


class Graphflow : public matching {
public:
    // a list of matching orders starting from each query edge
    // the first matching order also applies to the initial matching
    std::vector<std::vector<uint> > order_vs_; //????????
    std::vector<std::vector<uint> > order_csrs_;//???????????
    std::vector<std::vector<uint> > order_offs_;//??????????
    std::vector<std::vector<uint>> order_vertex_index;//???????????????§Ö?¦Ë?? order_vertex_index[0][u1]??????????????§µ?u1??????¦Ë??
//        std::unordered_map<std::pair<uint,uint>,std::vector<MatchRecord*>,pair_hash> edgeMaps;//????????????????????
    // uint s,t;//findMatch???????
    //???top k???????
    std::vector<MatchRecord *> topKSet;
    //?????????????
    std::vector<std::vector<StarGraph *>> globalStarIndex;//??????????????????????????????
    std::vector<SingleCandidate> match;//????????????  vertex density
    std::vector<std::vector<SingleCandidate>> matchCandidate;//id density
    std::vector<float> suffixMax;
    std::vector<float> isolatedMax;
    std::vector<std::vector<std::vector<Neighbor>>> rightNeighbor;//??????????id??
    //std::vector<LocalIndex>queryLocalIndexs;
    std::vector<std::vector<std::vector<int>>> globalVkMatchUk;//<vk,ak,uk>
    std::vector<std::vector<uint>> labelToQueryVertex;//???????????????
    std::vector<uint> queryVertexIndexInlabel;//??????????label?????§Ö???????
    std::vector<float> LocalStarIndex;//???????
    std::vector<std::vector<int>> matchLeftNeighborSum;//???§ß???????????
    std::vector<std::vector<int>> matchVetexLeftNeighbor;//???????????????????????
    std::vector<vector<float>> matchVetexSumweight;//??????????????????
    std::vector<std::vector<int>> leftNeighborIdSum;//???????????id??
    bool isUpdateIntopkset = false;
    int numAddTopk = 0;
    int allMatchFind = 0;
    long sumAllMatchFind = 0;
    long sumDeleteallMatchFind = 0;
    int numupdatestar = 0;
    long IsearchSpace=0,DsearchSpace=0,IdeterminCandite=0,DdeterminCandite=0;
    /* Timer total_search_time, total_print_time, total_densityFilter_time, total_update_globalIndex_time, total_updaterightNeighborCandidate_time,
             total_delete_time, total_delete_update_time;*/

public:
    Graphflow(Graph &query_graph, Graph &data_grasph, uint max_num_results,
              bool print_prep, bool print_enum, bool homo);

    ~Graphflow() override;

    void Preprocessing() override;

    void InitialMatching(const std::string &path) override;

    void AddEdge(uint v1, uint v2, uint label, float weight, uint timestamp) override;

    void AddEdgeWithGlobalIndex(uint v1, uint v2, uint label, float weight, uint timestamp) override;

    void RemoveEdge(uint v1, uint v2, uint label) override;

    void AddVertex(uint id, uint label) override;

    void RemoveVertex(uint id) override;

    void InitialTopK(const std::string &path) override;//????????????Top k???????

    void updateTopK() override;

    void deleteEdge(uint v1, uint v2) override;

    void deleteUpdateTopK() override;

    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

    void PrintAverageTime(int len);

private:
    void GenerateMatchingOrder();

    void FindMatches(uint flag, uint order_index, uint depth,
                     std::vector<uint> m, size_t &num_results, float density_s,
                     searchType type); //flag==0 initial flag==1 update
    int addMatchRecords(MatchRecord *r);

    void addStarGraph(StarGraph *s);

    void CreateStarIndex();

    float GetBackWeight(uint order_index, uint depth);

    void updateStarIndex(uint match_index, uint caddidate_v, const std::vector<uint> &canditeQueryVertexs);

    void updateStarIndex(uint match_index, uint caddidate_v, uint candidate_u, int candidate_v_index);

    void updateStarIndex(uint match_index, uint caddidate_v, uint candidate_u);

    vector<int> EdgeisInMatchOrder(Edge *edge);

    vector<int> EdgeisInMatchOrder(uint v1, uint v2, uint v1label, uint v2label, uint velabel);

    void searchMatches(int depth, uint matchorderindex, searchType flag);

    bool LabelFilter(uint data_v, uint query_v);

    void matchVertex(bool isFirstEdge, uint depth, uint data_v, float w);

    void matchVertex(int depth);

    void popVertex(uint depth, uint data_v);

    void popVertex(uint data_v, uint matchorderindex, uint depth, const std::vector<Neighbor> &uk_neighbor);

    void popVertex(uint depth);

    void densityFilter(uint matchorder_index, uint depth, std::vector<SingleCandidate> &singleVertexCandidate);

    void combination_helper(std::vector<std::vector<int>> &result, std::vector<int> &current,
                            const std::vector<std::vector<uint>> &nums, int k);

    std::vector<std::vector<int>> combination(const std::vector<std::vector<uint>> &nums);

    bool LDVertexCandidateCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> &needToIntersection,
                                std::vector<uint> &intersectresult);

    void setSingleVertexByIntersectionResult(std::vector<tuple<int, int, float>> &singleVertex,
                                             std::vector<uint> &intersectresult, std::vector<int> &r);

    void setIsolateVertexMatchResult(std::vector<int> &r, std::vector<int> &isolateVertex, float density);

    void setBatchVisited(std::vector<int> &r, bool flag);

    void recoverIsolateVertexMatchResult(std::vector<int> &IsolateVertexs);

    void sychronizeSingleVertexAndCandidate(std::vector<tuple<int, int, float>> &singleVertex,
                                            std::vector<uint> &intersectresult);

    void addMatchResult(uint matchorderindex, searchType type);

    std::vector<tuple<std::vector<int>, int, float>>
    combinationMatchResult(std::vector<std::vector<tuple<int, int, float>>> combinezIsolateVertexs);

    void combinationMatchResultHelp(std::vector<tuple<std::vector<int>, int, float>> &result, std::vector<int> &current,
                                    std::vector<std::vector<tuple<int, int, float>>> &combinezIsolateVertexs, int k,
                                    int tmin, float density
    );

    float findWeightBeforeIsolated();

    void CatesianProductWithIndex(int matchorderindex, searchType type, int curIndex, int depth, int len, int *hash,
                                  std::vector<std::vector<SingleCandidate>> &combinezIsolateVertexs,
                                  std::vector<int> &isolateVertexs, float &weight);

    int findTboundMaxIndex(float *Tbound, int *hash, int *nocan,
                           std::vector<std::vector<SingleCandidate>> &combinezIsolateVertexs, int len);

    bool isnoNextVertex(int *noscan, int len);

    void addMatchResultWithHeap(uint matchorderindex, searchType type);

    void CatesianProductWithHeap(int matchorderindex, searchType type, int depth, int len, int *hash,
                                 std::vector<std::vector<SingleCandidate>> &combinezIsolateVertexs,
                                 std::vector<int> &isolateVertexs, std::vector<int> &isolatedIndex, float &weight);

    void createLabelToQueryVertex();

    bool updaterightNeighborCandidate(int matchorderindex, uint uk, uint uk_neigh, bool isFirstEdge, uint vk,
                                      const std::vector<Neighbor> &uk_neighbor);

    void InitialLocalIndex(int matchorderindex);

    void getIntersetSingleCandidate(std::vector<SingleCandidate> &candidates, int matchorderindex, int depth);

    void deleteUpdateglobalVertexStarIndex(uint u1, uint v1, uint n);

    void deleteUpdateStarIndex(uint v1, uint v2, std::vector<int> &match);

    bool deleteMatchRecordWithEdge(uint v1, uint v1label, uint v2, uint v2label, uint label, std::vector<int> &match);

    void SearchMatchesWithEdge(uint m,uint v1,uint v2,uint weight,uint u1,uint u2,searchType type);

    void SearchMatchesWithGlobalIndexEdge(uint m, uint v1, uint v2, uint weight, uint u1, uint u2, searchType type);
    void searchMatchesWithGLobalIndex(int depth, uint matchorderindex, searchType flag) ;


};

#endif //MATCHING_GRAPHFLOW
