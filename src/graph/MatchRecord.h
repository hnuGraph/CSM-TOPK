//
// Created by ¸ßÉ­É­ on 2022/11/22.
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
    float density;
    std::vector<uint> vetexs=std::vector<uint>();//vertex set
public:
    MatchRecord(){};
    MatchRecord(float density_,std::vector<uint>vetexs_);
    ~MatchRecord(){};
    void AddVetex(uint u);
    void setDensity(float d);
    float getDensity();
    const std::vector<uint> &getVetex();
    std::string toString();
    std::string printMatchRecord();//print density tmin vetexs
    bool operator>(MatchRecord&m);
    bool operator==(const MatchRecord&m)const;
};


#endif //BASELINE_MATCHRECORD_H
