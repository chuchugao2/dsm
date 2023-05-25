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


/*class Vertex{
private:
    uint vertexId;
    uint tmin;
    float sumWeight;
public:
    Vertex(uint vertexId_,uint tmin_,uint sumWeight_):vertexId(vertexId_),tmin(tmin_),sumWeight(sumWeight_){};
    const bool operator>(const Vertex &v) const{
        if(this->sumWeight!=v.sumWeight){
            return this->sumWeight>v.sumWeight;
        }
        else if(this->vertexId!=v.vertexId){
            return this->vertexId>v.vertexId;
        }
    }
    const float GetSumWeight() const{
        return sumWeight;
    }
    const uint GetVertexId()const{
        return vertexId;
    }
};*/

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
    std::vector<MatchRecord*> allMatchRecords;
    std::vector<std::vector<StarGraph*>>qForwardNeighbors;//记录每种匹配下每个节点前向邻居以及最大权值
    std::vector<tuple<int,int,float>>match;//每个节点的匹配结果  vertex,tmin,density
    std::vector<std::vector<tuple<int,int,float>>>matchCandidate;//匹配序中，每个顶点的与其前邻居的密度和
    std::vector<std::vector<uint>>matchVertexCandidate;



public:
    Graphflow(Graph& query_graph, Graph& data_grasph, uint max_num_results,
              bool print_prep, bool print_enum, bool homo);
    ~Graphflow() override {};

    void Preprocessing() override;
    void InitialMatching(const std::string &path) override;

    void AddEdge(uint v1, uint v2, uint label,float weight,uint timestamp) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    void InitialTopK(const std::string &path) override;//得到初始化之后的Top k结果集合
    void updateTopK(uint num) override;
    void deleteEdge(uint v1,uint v2) override;
    void deleteUpdateTopK() override;
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;

private:
    void GenerateMatchingOrder();
    void FindMatches(uint flag,uint order_index, uint depth,
                     std::vector<uint> m, size_t &num_results,float density_s,uint tmin); //flag==0 initial flag==1 update
    bool addMatchRecords(MatchRecord* r);
    void addStarGraph(StarGraph *s);
    void CreateStarIndex();
    float GetBackWeight(uint order_index,uint depth);
    void updateStarIndex(uint match_index,uint order_vertex_index,uint candidate_u);
    vector<int> EdgeisInMatchOrder(Edge *edge);
    vector<int> EdgeisInMatchOrder(uint v1,uint v2,uint v1label,uint v2label,uint velabel);
    void searchMatches(uint matchorderindex,searchType flag);
    bool LabelFilter(uint data_v,uint query_v);
    void matchVertex(uint data_v,int tmin,float density);
    void matchVertex(uint data_v,int index,uint depth,std::vector<tuple<int,int,float>>&singleVertexCandidate);
    void matchVertex(std::vector<uint> & vertexCandidate,std::vector<tuple<int,int,float>>&singleVertexCandidate);
    void popVertex(uint data_v);
    void popVertex(uint data_v,int index,uint depth);
    void popVertex();
   // void densityFilter(uint matchorder_index,uint depth,std::vector<Vertex>&singleVertexCandidate,std::vector<uint>& candidate);
    void densityFilter(uint matchorder_index,uint depth,std::vector<uint>& candidate, std::vector<tuple<int,int,float>>&singleVertexCandidate,bool isSychronized);
    void combination_helper(std::vector<std::vector<int>>& result, std::vector<int>& current, const std::vector<std::vector<uint>>& nums, int k);
    std::vector<std::vector<int>> combination(const std::vector<std::vector<uint>>& nums);
    bool LDVertexCandidateCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection,  std::vector<uint> &intersectresult);
    void setSingleVertexByIntersectionResult( std::vector<tuple<int,int,float>>&singleVertex,std::vector<uint> &intersectresult,std::vector<int>&r);
    void setLDVertexMatchResult(std::vector<int>&r,std::vector<uint>&LDVertexs);
    void setIsolateVertexMatchResult(std::vector<int>&r,std::vector<int>&isolateVertex,uint tmin,float density);
    void setIsolateVertexMatchResult(std::vector<pair<int,int>>&r,std::vector<int>&isolateVertex,uint tmin,float density);
    void setBatchVisited(std::vector<int>&r,bool flag);
    void recoverLDVertexMatchResult(std::vector<int>&r,std::vector<uint>&LDVertexs);
    void recoverIsolateVertexMatchResult(std::vector<int>&IsolateVertexs);
    void sychronizeSingleVertexAndCandidate( std::vector<tuple<int,int,float>>&singleVertex,std::vector<uint> &intersectresult);
    void addMatchResult(uint matchorderindex,searchType type);
    std::vector<tuple<std::vector<int>,int,float>>combinationMatchResult(std::vector<std::vector<tuple<int,int,float>>>combinezIsolateVertexs);
    void combinationMatchResultHelp(std::vector<tuple<std::vector<int>,int,float>>&result,std::vector<int>&current,
                                    std::vector<std::vector<tuple<int,int,float>>>&combinezIsolateVertexs,int k,int tmin,float density
    );
    std::pair<int,float>findWeightAndTminBeforeIsolated();
    void CatesianProductWithIndex(int matchorderindex,searchType type,int curIndex,int depth,int len,int*hash,std::vector<std::vector<tuple<int,int,float>>>&combinezIsolateVertexs,std::vector<int>&isolateVertexs,int &tmin,float &weight);
    int findTboundMaxIndex(float *Tbound,int*hash,int*nocan,std::vector<std::vector<tuple<int,int,float>>>&combinezIsolateVertexs,int len);
    bool isnoNextVertex(int*noscan,int len);

};

#endif //MATCHING_GRAPHFLOW
