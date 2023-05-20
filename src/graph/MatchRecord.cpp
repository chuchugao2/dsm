//
// Created by ¸ßÉ­É­ on 2022/11/22.
//

#include "MatchRecord.h"


MatchRecord::MatchRecord(float density_, uint tmin_, std::vector <uint> vetexs_)
 :density(density_),tmin(tmin_),vetexs(vetexs_)
{}
float MatchRecord::getDensity() {
    return density;
}
void MatchRecord::setDensity(float d) {
    density=d;
}
void MatchRecord::setTmin(uint t) {
    tmin=t;
}
uint MatchRecord::getTmin() {
    return tmin;
}
std::vector<uint>* MatchRecord::getVetex() {
    return &vetexs;
}
void MatchRecord::AddVetex(uint u) {
    vetexs.emplace_back(u);
}

std::string MatchRecord::toString() {
    std::string str="density:";
    str+=std::to_string(density);
    str+=" tmin:";
    str+=std::to_string(tmin);
    str+=" vetexs:";

    for(uint i:vetexs){
        str+=std::to_string(i);
        str+=" ";
    }
    str+="\n";
    return str;
}
std::string MatchRecord::printMatchRecord() {
    std::string str="t ";
    str+=std::to_string(density);
    str+=" ";
    str+=std::to_string(tmin);
    str+=" ";
    for(int i=0;i<vetexs.size();i++){
        str+=std::to_string(vetexs[i]);
        if(i!=vetexs.size()-1){
            str+=" ";
        }
    }
    str+="\n";
    return str;
}
bool MatchRecord::operator>(MatchRecord &m) {
    if(this->density!=m.density){
        return this->density>m.density;
    }else if(this->tmin!=m.tmin){
        return this->tmin>m.tmin;
    }else {
        return this->vetexs>m.vetexs;
    }
}
bool MatchRecord::operator==(MatchRecord &m) {
    return density==m.density&&tmin==m.tmin&&vetexs==m.vetexs;
}

