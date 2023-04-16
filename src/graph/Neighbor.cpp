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
    else{
        return this->edgeWeight>m.edgeWeight;
    }
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