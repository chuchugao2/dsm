//
// Created by �߳��� on 2023/3/23.
//

#include "StarGraph.h"
void StarGraph::SetStarMaxWeight(uint w) {
    maxWeight=w;
}
const uint StarGraph::GetStarMaxWeight() {
    return maxWeight;
}
void StarGraph::InitalmaxWeight() {
    maxWeight=queryVertex.size()*mw;
}
const uint StarGraph::GetForwardNeighborNum() {
    return queryVertex.size();
}
StarGraph::~StarGraph() {
    for(int i=0;i<queryVertex.size();i++){
        delete queryVertex[i];
        queryVertex[i]=NULL;
    }
    this->queryVertex.clear();
}
void StarGraph::computeMaxWeight() {
    maxWeight=0;
    for(auto f:queryVertex){
        maxWeight+=f->GetMaxWeight();
    }
}

std::vector<ForwardNeighbor *> &StarGraph::GetqueryVertex() {
return queryVertex;
}