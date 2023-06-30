//
// Created by �߳��� on 2023/6/18.
//

#include "Subgraph.h"



Subgraph::Subgraph(uint vertexSize,uint queryVertexSize):
vNeighbors(vertexSize),
matchCandidate(queryVertexSize),
setVNeighbors(vertexSize)
{}

void Subgraph::addQueryVertexCandidate(uint q, uint v) {
    //�������
    auto it = std::lower_bound(matchCandidate[q].begin(), matchCandidate[q].end(), v);
    // ��ȷ����λ�ò����µ�Ԫ��
    matchCandidate[q].insert(it, v);
//    matchCandidate[q].emplace_back(v);
}
bool Subgraph::AddEdge(uint u1,uint u2,uint v1, uint v1label,uint v2,uint v2label, uint label, float weight) {
//�ҵ����ڵ���v2�ĵ�������λ��
    Neighbor neighbor(v2,v2label,label,weight,u2,u1);
    if(setVNeighbors[v1].count(neighbor)){
        return false;
    }
    else{
        vNeighbors[v1].emplace_back(neighbor);
        setVNeighbors[v1].insert(neighbor);
    }
    Neighbor neighbor2(v1, v1label,label,weight,u1,u2);
    if(setVNeighbors[v2].count(neighbor2)){
        return false;
    }
    else{
        vNeighbors[v2].emplace_back(neighbor2);
        setVNeighbors[v2].insert(neighbor2);
    }
    edge_count_++;
    return true;
}

bool Subgraph::RemoveEdge(uint v1, uint v2,uint u1,uint u2)
{
    Neighbor neighbor(v2,u2,u1);
    if(!setVNeighbors[v1].count(neighbor)){
        return false;
    }
    for(auto it=vNeighbors[v1].begin();it!=vNeighbors[v1].end();it++){
        if((*it)==neighbor){
            vNeighbors[v1].erase(it);
            break;
        }
    }
    setVNeighbors[v1].erase(neighbor);

    Neighbor neighbor2(v1,u1,u2);
    if(!setVNeighbors[v2].count(neighbor2)){
        return false;
    }
    for(auto it=vNeighbors[v2].begin();it!=vNeighbors[v2].end();it++){
        if((*it)==neighbor2){
            vNeighbors[v2].erase(it);
            break;
        }
    }
    setVNeighbors[v2].erase(neighbor2);
    edge_count_--;
    return true;
}
void Subgraph::deleteQueryVertexCandidate(uint q, uint v) {
    std::vector<uint>&candidate=matchCandidate[q];
    auto it=std::lower_bound(candidate.begin(),candidate.end(),v);

    if (it != candidate.end() && *it ==v) {
        // ɾ��Ŀ��Ԫ��
        candidate.erase(it);
    }


}
const std::vector<Neighbor>& Subgraph::GetVNeighbors(uint v) const {
    return vNeighbors[v];
}



