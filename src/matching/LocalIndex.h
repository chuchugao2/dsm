//
// Created by ¸ß³þ³þ on 2023/5/31.
//

#ifndef BASELINE_LOCALINDEX_H
#define BASELINE_LOCALINDEX_H

#include <vector>
#include <cfloat>
#include "limits.h"

class LocalIndex {
protected:
    std::vector<std::pair<int,float>>candidate;//data vertex/weight
    float maxWeight=FLT_MAX;
    int maxId=INT_MAX;
    bool isFirst= true;
public:
    LocalIndex(){};
    ~LocalIndex(){};
    std::vector<std::pair<int,float>>& getCandidate();
    void setMaxWeght(float max_weight);
    void setmaxId(int max_id);
    float getMaxWeight();
    int getMaxId();
    void insertCandidate(int id,float weight);
    void clearCandidate();
    void setIsFirst(bool flag);
    bool getIsFirst();
};

#endif //BASELINE_LOCALINDEX_H
