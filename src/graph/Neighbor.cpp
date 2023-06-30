//
// Created by ¸ß³þ³þ on 2023/3/23.
//

#include "Neighbor.h"
uint Neighbor::getVertexId() const {
    return toVertexId;
}
 bool Neighbor::operator>(const Neighbor &m) const{
    if (this->edgeLabel!=m.edgeLabel)
    {
        return this->edgeLabel>m.edgeLabel;
    }
    else if(this->toVertexLabel!=m.toVertexLabel)
    {
        return this->toVertexLabel>m.toVertexLabel;
    }
    else if(this->fromVertexId!=m.fromVertexId){
        return this->fromVertexId>m.fromVertexId;
    }
    else if(this->matchQueryVertexId!=m.matchQueryVertexId){
        return this->matchQueryVertexId>m.matchQueryVertexId;
    }
    else{
        return this->edgeWeight>m.edgeWeight;
    }
}
bool Neighbor::operator==(const Neighbor &m) const {
        return toVertexId == m.toVertexId &&
               matchQueryVertexId == m.matchQueryVertexId &&
               fromVertexId == m.fromVertexId;
}

std::pair<uint,uint> Neighbor::GetelabelAndVertexLabel() const {
    return std::make_pair(edgeLabel,toVertexLabel);
}
float Neighbor::GetEdgeWeight() const {
    return edgeWeight;
}
uint Neighbor::getVertexLabel() const {
    return toVertexLabel;
}
uint Neighbor::GetEdgelabel() const {
    return edgeLabel;
}
const uint Neighbor::getMatchQueryVertexId() const{
    return matchQueryVertexId;
}
const uint Neighbor::getfromVertexId() const {
    return fromVertexId;
}