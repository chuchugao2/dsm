//
// Created by ¸ß³þ³þ on 2023/3/23.
//

#ifndef BASELINE_FORWARDNEIGHBOR_H
#define BASELINE_FORWARDNEIGHBOR_H
#include <vector>
#include "../utils/types.h"
#include "../utils/globals.h"

class ForwardNeighbor{
protected:
    uint toVertexIndex=0;
    uint toVertexId=0;
    uint toVertexLabel;
    uint edgeLabel;
    float maxWeight=mw;
public:
    ForwardNeighbor(){};
    ForwardNeighbor(uint toVertexIndex_,uint toVertexId_,uint toVertexLabel_,uint edgeLabel_,float maxWeight_):toVertexIndex(toVertexIndex_),toVertexId(toVertexId_),toVertexLabel(toVertexLabel_),edgeLabel(edgeLabel_),maxWeight(maxWeight_){};
    ForwardNeighbor(uint toVertexLabel_,uint edgeLabel_,float maxWeight_):toVertexLabel(toVertexLabel_),edgeLabel(edgeLabel_),maxWeight(maxWeight_){};
    ForwardNeighbor(uint toVetexId_,uint toVertexLabel_,uint edgeLabel_):toVertexId(toVetexId_),toVertexLabel(toVertexLabel_),edgeLabel(edgeLabel_){};
    ForwardNeighbor(uint toVertexIndex_,uint toVetexId_,uint toVertexLabel_,uint edgeLabel_):toVertexIndex(toVertexIndex_),toVertexId(toVetexId_),toVertexLabel(toVertexLabel_),edgeLabel(edgeLabel_){};
    ~ForwardNeighbor(){};
    const bool operator>(const ForwardNeighbor &f){
      if(this->edgeLabel!=f.edgeLabel)
          return this->edgeLabel>f.edgeLabel;
      else{
          return this->toVertexLabel>f.toVertexLabel;
      }
    };
    const float GetMaxWeight();
    const std::pair<uint,uint> GetelabelAndVertexLabel() ;
    const uint GetVetexId()const;
    const uint GetElabel()const;
    const uint GetVetexIndex()const;
};

#endif //BASELINE_FORWARDNEIGHBOR_H
