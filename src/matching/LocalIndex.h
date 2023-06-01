//
// Created by ¸ß³þ³þ on 2023/5/31.
//

#ifndef BASELINE_LOCALINDEX_H
#define BASELINE_LOCALINDEX_H

#include <vector>
#include <cfloat>
#include "limits.h"

class LocalIndex {
public:
    std::vector<std::pair<int,float>>candidate;//data vertex/weight
    float maxWeight=FLT_MAX;
    int maxIndex=INT_MAX;
public:
    LocalIndex(){};
    ~LocalIndex(){};
    std::vector<std::pair<int,float>>& getCandidate();
    void setMaxWeght(float max_weight);
    void setmaxIndex(int max_index);
    float getMaxWeight();
    int getMaxIndex();
    void insertCandidate(int id,float weight);
    void clearCandidate();
};

#endif //BASELINE_LOCALINDEX_H
