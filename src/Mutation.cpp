//
// Created by Smith, Maureen on 31.05.18.
//

#include "Mutation.hpp"

Mutation::Mutation(const unsigned int p, const unsigned int s):position(p), symbol(s){}

const unsigned int Mutation::getPosition() const {
    return position;
}

const unsigned int Mutation::getSymbol() const {
    return symbol;
}

//unsigned int Mutation::getCount() const {
//    return count;
//}
//
//double Mutation::getFreq() const {
//    return freq;
//}
//
//void Mutation::setCount(unsigned int count) {
//    Mutation::count = count;
//
//}
//
//void Mutation::setFreq(double freq) {
//    Mutation::freq = freq;
//}
//
//double Mutation::getNoise() const {
//    return noise;
//}
//
//void Mutation::setNoise(double noise) {
//    Mutation::noise = noise;
//}

bool Mutation::operator<(const Mutation &mut) const {
    return this->getPosition() < mut.getPosition();
}