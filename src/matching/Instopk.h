//
// Created by 高楚楚 on 2023/5/16.
//

#ifndef BASELINE_INSTOPK_H
#define BASELINE_INSTOPK_H


#include <map>
#include "matching.h"
#include "math.h"
class Path{
public:
    std::vector<int>nodes;
public:
    Path(std::vector<int>nodes_):nodes(nodes_){};
    Path(){};
    bool operator <(const Path & p)const{
        return nodes<p.nodes;
    }
};
class Instopk: public matching{
public:
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
    std::vector<std::vector<int>>matchVertexCandidate;//每个节点的候选节点
    std::vector<std::vector<std::map<std::string,int>>>TopologyIndex;//拓扑索引
    std::vector<std::vector<std::map<std::string,int>>>queryTopologyIndex;//查询图拓扑索引
    std::vector<std::vector<std::map<std::string,float>>>MNW;//最大源路径权重<d,id<<a,1><b,1><c,1>>
    std::map<std::string,std::vector<Edge>>sortEdgeList;//sortlist
    std::vector<std::vector<std::string>>dToOrderType;
    std::map<int,std::vector<std::string>>querydToOrderType;
    uint dist;
    std::map<std::string,std::vector<std::vector<int>>>node2EdgeListPointers;//key 0#0边类型，map[0]为v1在类型为0#0的索引个数
    std::map<std::string,int>pointers;//sortEdgelist指针
    std::set<std::string>queryEdgeTypes;//查询图的所有边类型
    std::map<std::string,std::vector<Edge>>queryEdgeType2Edges;//查询边每种边的类型对应的查询边id
    std::map<std::string,std::string>queryEdgesToEdgeType;//查询边id对应的查询边类型
    std::map<std::string,int>queryEdgeToIndex;//每条边对应的索引


public:
    Instopk(Graph& query_graph, Graph& data_grasph, uint max_num_results,
    bool print_prep, bool print_enum, bool homo,uint dist);
    ~Instopk() override {};
    void Preprocessing() override;
    void InitialMatching(const std::string &path) override;
    void AddEdge(uint v1, uint v2, uint label,float weight,uint timestamp) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    void InitialTopK(const std::string &path) override;//得到初始化之后的Top k结果集合
    void updateTopK() override;
    void deleteEdge(uint v1,uint v2) override;
    void deleteUpdateTopK() override;
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
private:
    void CreateSortEdgeList();//创建排序边列表
    void GenerateMatchingOrder();
    void CreateTopologyAndMNWIndex();
    void CreateQueryTopology();
    void SearchMatchesWithSortedList();//按照边搜索过程
    void InitialPointers();//初始化排序边列表指针
    std::pair<int,int>splitString(std::string s);
    float GetUperBoundWithPath(const std::vector<int>&consideredEdgeIndex,const std::vector<std::string>&pc);
    std::set<Path>GetPaths(int i,const std::vector<std::pair<int,int>>&coverEdges);
    bool addMatchRecords(MatchRecord* r);
    void InitialqueryCandidate();//初始化查询节点的候选节点
    bool isContainMatchRecord(MatchRecord *m);
    bool count(std::vector<int>&t,int e);
    bool count(  std::vector<Edge>&t,Edge e);
    void updateSortEdgelist(uint v1,uint v2);
    void updateMNWIndexAndDataTopologyIndex();
    void updateQueryCandidate();
    void insertEdgeToSortEdgelist(std::string &str,Edge &edge);



};


#endif //BASELINE_INSTOPK_H
