//
// Created by �߳��� on 2023/3/23.
//

#ifndef BASELINE_STARGRAPH_H
#define BASELINE_STARGRAPH_H

#include <vector>
#include "../utils/types.h"
#include "../utils/globals.h"
#include "ForwardNeighbor.h"


class StarGraph {
protected:
    std::vector<ForwardNeighbor *> queryVertex;
    float maxWeight;
    uint MatchDataVertexId = UINT_MAX;
public:
    StarGraph() {};

    StarGraph(std::vector<ForwardNeighbor *> q) {
        queryVertex = q;
    }

    ~StarGraph();

    void AddForwardNeighbor(ForwardNeighbor *f) {
        queryVertex.emplace_back(f);
    }

    void InitalmaxWeight();

    void computeMaxWeight();

    const uint getStarMaxWeight();

    const uint getMatchDataVertexId();

    void setMatchDataVertexId(uint id);

    void setStarMaxWeight(uint w);

    const uint GetForwardNeighborNum();

    std::vector<ForwardNeighbor *> &GetqueryVertex();
};

#endif //BASELINE_STARGRAPH_H
