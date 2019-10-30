//
// Created by Smith, Maureen on 31.05.18.
//

#include "Mutation.hpp"

Mutation::Mutation(const unsigned int p, const double kd):position(p), kd(kd) {}

const unsigned int Mutation::getPosition() const {
    return position;
}

const double Mutation::getKd() const {
    return kd;
}

unsigned int Mutation::getCount() const {
    return count;
}

double Mutation::getFreq() const {
    return freq;
}

void Mutation::setCount(unsigned int count) {
    Mutation::count = count;
}

void Mutation::setFreq(double freq) {
    Mutation::freq = freq;
}

double Mutation::getNoise() const {
    return noise;
}

void Mutation::setNoise(double noise) {
    Mutation::noise = noise;
};