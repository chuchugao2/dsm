//
// Created by ¸ß³þ³þ on 2023/7/20.
//

#ifndef BASELINE_SUBNEIGHBOR_H
#define BASELINE_SUBNEIGHBOR_H

#include <utility>
#include "../utils/globals.h"

class SubNeighbor {
    protected:
            uint toVertexId;
    uint toVertexLabel;
    uint edgeLabel;
    float edgeWeight;
    uint matchQueryVertexId;
    uint fromQueryVertexId;
public:
    SubNeighbor();
    SubNeighbor(uint toVertexId_,uint toVertexLabel_,uint edgeLabel_,float edgeWeight_,uint matchQueryVertexId_,uint fromVertexId_):
    toVertexId(toVertexId_),
    toVertexLabel(toVertexLabel_),
    edgeLabel(edgeLabel_),
    edgeWeight(edgeWeight_),
    matchQueryVertexId(matchQueryVertexId_),
    fromQueryVertexId(fromVertexId_)
    {};
    SubNeighbor(uint toVertexId_,uint toVertexLabel_,float edgeLabel_):
    toVertexId(toVertexId_),
            toVertexLabel(toVertexLabel_),
            edgeLabel(edgeLabel_),
            edgeWeight(0){};
    SubNeighbor(uint toVertexId_,uint matchQueryVertexId_,uint fromVertexId_):
    toVertexId(toVertexId_),
    matchQueryVertexId(matchQueryVertexId_),
    fromQueryVertexId(fromVertexId_){};
    uint getVertexId() const;
    std::pair<uint,uint> GetelabelAndVertexLabel() const;
    bool operator>(const SubNeighbor &m) const;
    bool operator!=(const SubNeighbor&m)const;
    bool operator==(const SubNeighbor&m)const;
    float GetEdgeWeight()const;
    uint getVertexLabel()const;
    uint GetEdgelabel()const;
    const uint getMatchQueryVertexId()const ;
    const uint getfromVertexId()const ;

};


#endif //BASELINE_SUBNEIGHBOR_H
