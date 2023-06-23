//
// Created by 高楚楚 on 2023/6/18.
//

#include "Subgraph.h"



Subgraph::Subgraph(uint vertexSize,uint queryVertexSize):
vNeighbors(vertexSize),
matchCandidate(queryVertexSize),
setVNeighbors(vertexSize)
{}

void Subgraph::addQueryVertexCandidate(uint q, uint v) {
    //按序插入
    auto it = std::lower_bound(matchCandidate[q].begin(), matchCandidate[q].end(), v);
    // 在确定的位置插入新的元素
    matchCandidate[q].insert(it, v);
//    matchCandidate[q].emplace_back(v);
}
bool Subgraph::AddEdge(uint u1,uint u2,uint v1, uint v1label,uint v2,uint v2label, uint label, float weight) {
//找到大于等于v2的迭代器的位置
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

void Subgraph::RemoveEdge(uint v1, uint v2,uint u1,uint u2)
{
    for(auto it=vNeighbors[v1].begin();it!=vNeighbors[v1].end();it++){
        if(it->getVertexId()==v2&&it->getMatchQueryVertexId()==u2&&it->getfromVertexId()==u1){
            vNeighbors[v1].erase(it);
        }
    }

    for(auto it=vNeighbors[v2].begin();it!=vNeighbors[v2].end();it++){
        if(it->getVertexId()==v1&&it->getMatchQueryVertexId()==u1&&it->getfromVertexId()==u2){
            vNeighbors[v2].erase(it);
        }
    }
    edge_count_--;
}
void Subgraph::deleteQueryVertexCandidate(uint q, uint v) {
    std::vector<uint>&candidate=matchCandidate[q];
    for(int i=0;i<candidate.size();i++){
        if(candidate[i]==v){
            candidate.erase(candidate.begin()+i);
            break;
        }
    }
}
const std::vector<Neighbor>& Subgraph::GetVNeighbors(uint v) const {
    return vNeighbors[v];
}



