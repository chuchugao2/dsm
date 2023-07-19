//
// Created by ¸ß³þ³þ on 2023/3/23.
//

#ifndef BASELINE_NEIGHBOR_H
#define BASELINE_NEIGHBOR_H

#include "../utils/globals.h"
#include <utility>

class Neighbor {
protected:
    uint toVertexId;
    uint toVertexLabel;
    uint edgeLabel;
    float edgeWeight;
    uint timetamp;
public:
    Neighbor();

    Neighbor(uint toVertexId_, uint toVertexLabel_, uint edgeLabel_, float edgeWeight_, uint timestap_) {
        this->toVertexId = toVertexId_;
        this->toVertexLabel = toVertexLabel_;
        this->edgeLabel = edgeLabel_;
        this->edgeWeight = edgeWeight_;
        this->timetamp = timestap_;
    };

    Neighbor(uint toVertexId_, uint toVertexLabel_, uint edgeLabel_) {
        this->toVertexId = toVertexId_;
        this->toVertexLabel = toVertexLabel_;
        this->edgeLabel = edgeLabel_;
        this->edgeWeight = 0;
        this->timetamp = 0;
    };
    Neighbor(uint toVertexLabel_, uint edgeLabel_) {
        this->toVertexLabel = toVertexLabel_;
        this->edgeLabel = edgeLabel_;
    };

    Neighbor(uint toVertexId_, uint toVertexLabel_, uint edgeLabel_, float edgeWeight_) {
        this->toVertexId = toVertexId_;
        this->toVertexLabel = toVertexLabel_;
        this->edgeLabel = edgeLabel_;
        this->edgeWeight = edgeWeight_;
    };

    uint getVertexId() const;

    std::pair<uint, uint> GetelabelAndVertexLabel() const;

    bool operator>(const Neighbor &m) const;

    bool operator==(const Neighbor &m) const;

    bool operator!=(const Neighbor &m) const;

    float GetEdgeWeight() const;

    uint getVertexLabel() const;

    uint GetEdgelabel() const;

};

#endif //BASELINE_NEIGHBOR_H
