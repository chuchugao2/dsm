//
// Created by �߳��� on 2023/3/23.
//

#include "ForwardNeighbor.h"

const std::pair<uint, uint> ForwardNeighbor::GetelabelAndVertexLabel() {
    return std::make_pair(edgeLabel, toVertexLabel);
}

const float ForwardNeighbor::GetMaxWeight() {
    return maxWeight;
}

const uint ForwardNeighbor::GetVetexId() const {
    return toVertexId;
}

const uint ForwardNeighbor::GetElabel() const {
    return edgeLabel;
}
const uint ForwardNeighbor::GetVertexLabel() const {
    return toVertexLabel;
}
void ForwardNeighbor::setMaxWeight(float m) {
    maxWeight=m;
}
void ForwardNeighbor::setMatchDataVertexId(uint id) {
    MatchDataVertexId=id;
}
const uint ForwardNeighbor::getMatchDataVertexId() {
    return MatchDataVertexId;
}

