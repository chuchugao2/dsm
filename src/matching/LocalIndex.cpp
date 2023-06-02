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
void LocalIndex::setmaxId(int max_id) {
    maxId=max_id;
}
float LocalIndex::getMaxWeight() {
    return maxWeight;
}
int LocalIndex::getMaxId() {
    return maxId;
}
void LocalIndex::insertCandidate(int id, float weight) {
    candidate.emplace_back(id,weight);
}
void LocalIndex::clearCandidate() {
    candidate.resize(0);
    maxWeight=FLT_MAX;
    maxId=INT_MAX;
    isFirst= true;
}
void LocalIndex::setIsFirst(bool flag) {
    isFirst= flag;
}
bool LocalIndex::getIsFirst() {
    return isFirst;
}