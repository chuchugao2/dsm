//
// Created by ��ɭɭ on 2022/11/22.
//

#ifndef BASELINE_MATCHRECORD_H
#define BASELINE_MATCHRECORD_H
#include "sys/types.h"
#include "iostream"
#include "vector"
#include "../utils/globals.h"
#include "algorithm"
class MatchRecord {
protected:
    float density;//�ܶ�
    uint tmin;//��С��ʱ���
    std::vector<uint> vetexs=std::vector<uint>();//�ڵ㼯��
public:
    MatchRecord(){};
    MatchRecord(float density_, uint tmin_, std::vector<uint>vetexs_);
    ~MatchRecord(){};
    void AddVetex(uint u);
    void setDensity(float d);
    void setTmin(uint t);
    uint getTmin();
    float getDensity();
    std::vector<uint> *getVetex();
    std::string toString();
    std::string printMatchRecord();//print density tmin vetexs
    bool operator>(MatchRecord&m);
    bool operator==(MatchRecord&m);
};


#endif //BASELINE_MATCHRECORD_H
