//
// Created by ¸ß³þ³þ on 2023/7/20.
//

#include "SubNeighbor.h"
uint SubNeighbor::getVertexId() const {
    return toVertexId;
}
bool SubNeighbor::operator>(const SubNeighbor &m) const{
    if (this->edgeLabel!=m.edgeLabel)
    {
        return this->edgeLabel>m.edgeLabel;
    }
    else if(this->toVertexLabel!=m.toVertexLabel)
    {
        return this->toVertexLabel>m.toVertexLabel;
    }
    else if(this->fromQueryVertexId!=m.fromQueryVertexId){
        return this->fromQueryVertexId>m.fromQueryVertexId;
    }
    else if(this->matchQueryVertexId!=m.matchQueryVertexId){
        return this->matchQueryVertexId>m.matchQueryVertexId;
    }
    else if((int)edgeWeight!=(int)m.edgeWeight){
        return this->edgeWeight>m.edgeWeight;
    }
    else{
        return this->toVertexId>m.toVertexId;
    }
}
bool SubNeighbor::operator==(const SubNeighbor &m) const {
    return toVertexId == m.toVertexId &&
           matchQueryVertexId == m.matchQueryVertexId &&
           fromQueryVertexId == m.fromQueryVertexId;
}
bool SubNeighbor::operator!=(const SubNeighbor &m) const {
    return !(*this==m);
}

std::pair<uint,uint> SubNeighbor::GetelabelAndVertexLabel() const {
    return std::make_pair(edgeLabel,toVertexLabel);
}
float SubNeighbor::GetEdgeWeight() const {
    return edgeWeight;
}
uint SubNeighbor::getVertexLabel() const {
    return toVertexLabel;
}
uint SubNeighbor::GetEdgelabel() const {
    return edgeLabel;
}
const uint SubNeighbor::getMatchQueryVertexId() const{
    return matchQueryVertexId;
}
const uint SubNeighbor::getfromVertexId() const {
    return fromQueryVertexId;
}