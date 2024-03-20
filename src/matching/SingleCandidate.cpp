//
// Created by ¸ß³þ³þ on 2023/6/3.
//

#include "SingleCandidate.h"
void SingleCandidate::setVertexId(int vertexID_) {
    vertexId=vertexID_;
}
void SingleCandidate::setSumWeight(float sumWeight_) {
    sumWeight=sumWeight_;
}
const int &SingleCandidate::getVertexId() const {
    return vertexId;
}

const float &SingleCandidate::getSumWeight() const {
    return sumWeight;
}
void SingleCandidate::clearSingleCandidate() {
    vertexId=0;
    sumWeight=0;
}
void SingleCandidate::setIsolateSingleCandate() {
    vertexId=-1;
    sumWeight=-1;
}
bool SingleCandidate::operator<(const SingleCandidate &s) const {
    if(sumWeight!=s.sumWeight)
        return sumWeight>s.sumWeight;
    else{
        return vertexId>s.vertexId;
    }
}
void SingleCandidate::setFlag(int i) {
    flag=i;
}
const int SingleCandidate::getFlag() const {
    return flag;
}
void SingleCandidate::addSumWeight(float weight) {
    sumWeight+=weight;
}
void SingleCandidate::addFlag() {
    flag++;
}