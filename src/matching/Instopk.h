//
// Created by �߳��� on 2023/5/16.
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
    std::vector<std::vector<uint> > order_vs_; //ƥ��ĵ�˳��
    std::vector<std::vector<uint> > order_csrs_;//ƥ��ڵ�ǰ���ھ�
    std::vector<std::vector<uint> > order_offs_;//ƥ��ڵ������
    std::vector<std::vector<uint>> order_vertex_index;//ÿ���ڵ��ڸ�ƥ�����е�λ�� order_vertex_index[0][u1]��ʾ��һ��ƥ�����У�u1������λ��
//        std::unordered_map<std::pair<uint,uint>,std::vector<MatchRecord*>,pair_hash> edgeMaps;//��Ӧ��ÿ���ߵĸ��������ṹ
    // uint s,t;//findMatchʱ���¼�ֵ
    //��¼top k��¼����
    std::vector<MatchRecord*> topKSet;
    //��¼����ƥ��Ľ��
    std::vector<MatchRecord*> allMatchRecords;
    std::vector<std::vector<StarGraph*>>qForwardNeighbors;//��¼ÿ��ƥ����ÿ���ڵ�ǰ���ھ��Լ����Ȩֵ
    std::vector<std::tuple<int,int,float>>match;//ÿ���ڵ��ƥ����  vertex,tmin,density
    std::vector<std::vector<std::tuple<int,int,float>>>matchCandidate;//ƥ�����У�ÿ�����������ǰ�ھӵ��ܶȺ�
    std::vector<std::vector<int>>matchVertexCandidate;//ÿ���ڵ�ĺ�ѡ�ڵ�
    std::map<int,std::map<int,std::map<std::string,int>>>TopologyIndex;//��������
    std::map<int,std::map<int,std::map<std::string,int>>>queryTopologyIndex;//��ѯͼ��������
    std::map<int,std::map<int,std::map<std::string,float>>>MNW;//���Դ·��Ȩ��<d,id<<a,1><b,1><c,1>>
    std::map<std::string,std::vector<Edge>>sortEdgeList;//sortlist
    std::map<int,std::vector<std::string>>dToOrderType;
    std::map<int,std::vector<std::string>>querydToOrderType;
    uint dist;
    std::map<std::string,std::map<int,std::vector<int>>>node2EdgeListPointers;
    std::map<std::string,int>pointers;//sortEdgelistָ��
    std::set<std::string>queryEdgeTypes;//��ѯͼ�����б�����
    std::map<std::string,std::vector<Edge>>queryEdgeType2Edges;//��ѯ��ÿ�ֱߵ����Ͷ�Ӧ�Ĳ�ѯ��id
    std::map<std::string,std::string>queryEdgesToEdgeType;//��ѯ��id��Ӧ�Ĳ�ѯ������
    std::map<std::string,int>queryEdgeToIndex;//ÿ���߶�Ӧ������


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
    void InitialTopK(const std::string &path) override;//�õ���ʼ��֮���Top k�������
    void updateTopK(uint num) override;
    void deleteEdge(uint v1,uint v2) override;
    void deleteUpdateTopK() override;
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
private:
    void CreateSortEdgeList();//����������б�
    void GenerateMatchingOrder();
    void CreateTopologyAndMNWIndex();
    void CreateQueryTopology();
    void SearchMatchesWithSortedList();//���ձ���������
    void InitialPointers();//��ʼ��������б�ָ��
    std::pair<int,int>splitString(std::string s);
    float GetUperBoundWithPath(std::set<int>consideredEdgeIndex,std::vector<std::string>pc);
    std::set<Path>GetPaths(int i,std::vector<std::pair<int,int>>coverEdges);
    bool addMatchRecords(MatchRecord* r);
    void InitialqueryCandidate();//��ʼ����ѯ�ڵ�ĺ�ѡ�ڵ�
    bool isContainMatchRecord(MatchRecord *m);

};


#endif //BASELINE_INSTOPK_H
