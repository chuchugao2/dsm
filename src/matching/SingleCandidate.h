//
// Created by ¸ß³þ³þ on 2023/6/3.
//

#ifndef BASELINE_SINGLECANDIDATE_H
#define BASELINE_SINGLECANDIDATE_H


class SingleCandidate {
protected:
    int vertexId=0;
    float sumWeight=0;
public:
    SingleCandidate(){};
    SingleCandidate(int vertexID_,float sumWeight_):vertexId(vertexID_),sumWeight(sumWeight_){};
    void setVertexId(int vertexID_);
    void setSumWeight(float sumWeight_);
    const int &getVertexId()const;
    const float &getSumWeight()const;
    ~SingleCandidate(){};
    void clearSingleCandidate();
    void setIsolateSingleCandate();
    bool operator<(const SingleCandidate&s)const;
};


#endif //BASELINE_SINGLECANDIDATE_H
