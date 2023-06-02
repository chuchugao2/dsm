//
// Created by ¸ß³þ³þ on 2023/3/23.
//

#include "ForwardNeighbor.h"
const std::pair<uint,uint> ForwardNeighbor::GetelabelAndVertexLabel()  {
    return std::make_pair(edgeLabel,toVertexLabel);
}
const float ForwardNeighbor::GetMaxWeight() {
    return maxWeight;
}
 const uint ForwardNeighbor::GetVetexId() {
    return toVetexId;
}
const uint ForwardNeighbor::GetElabel(){
    return edgeLabel;
}