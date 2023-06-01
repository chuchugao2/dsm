//
// Created by ¸ß³þ³þ on 2023/5/31.
//

#include "LocalIndex.h"
std::vector<std::pair<int,float>>& LocalIndex::getCandidate() {
    return candidate;
}
void LocalIndex::setMaxWeght(float max_weight) {
    maxWeight=max_weight;
}
void LocalIndex::setmaxIndex(int max_index) {
    maxIndex=max_index;
}
float LocalIndex::getMaxWeight() {
    return maxWeight;
}
int LocalIndex::getMaxIndex() {
    return maxIndex;
}
void LocalIndex::insertCandidate(int id, float weight) {
    candidate.push_back(std::make_pair(id,weight));
}
void LocalIndex::clearCandidate() {
    candidate.resize(0);
    maxWeight=FLT_MAX;
    maxIndex=INT_MAX;
}