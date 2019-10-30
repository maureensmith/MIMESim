//
// Created by Smith, Maureen on 31.05.18.
//

#include "FunctionalSequence.hpp"
#include "Constants.hpp"
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>

std::vector<double> FunctionalSequence::drawKdValues() {
    auto& constants = constants::Constants::get_instance();
    std::vector<double> kds(constants.L);
    const auto seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
    std::default_random_engine generator(seed);
    std::bernoulli_distribution bd(constants.P_EFFECT);
    std::lognormal_distribution<double> lnd(0,1);
    for(int i=0; i<kds.size(); ++i) {
        //first sample the positions which have an effect (binomial/bernoulli distributed), then the strength of the effect(log normal distributed)
        kds[i] = bd(generator)?lnd(generator):constants.KD_WT;
    }
    return kds;
}

std::vector<double> FunctionalSequence::drawEpistasis() {
    auto& constants = constants::Constants::get_instance();
    std::vector<double> epistasis(constants.PW);
    const auto seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
    std::default_random_engine generator(seed);
    std::bernoulli_distribution bd(constants.P_EPISTASIS);
    std::lognormal_distribution<double> lnd(0,1);
    //TODO: Kann epistasis auch vorhanden sein, wenn eine einzelmutation keinen effect hat? (hatte es in R so realisiert dass beide einene Effekt haben müssen.
    for(int i=0; i<epistasis.size(); ++i) {
        // first sample if position pair has epistatic effect (bernoulli) and then the value of the epistasis (log normal distributed)
        epistasis[i] = bd(generator)?lnd(generator):constants.NO_EPISTASIS;
    }
    return epistasis;
}

FunctionalSequence& FunctionalSequence::get_instance(){
    static FunctionalSequence instance;
    return instance;
}

const double &FunctionalSequence::getKd(const unsigned int i) {
    //TODO abfangen
    if(i>0 && i<=kds.size())
        return kds[i-1];
    else
       std::cerr << "TODO: Fehler Abfangen";
}

unsigned int FunctionalSequence::getMatrixVectorIndex(unsigned int i, unsigned int j) {
    auto& c = constants::Constants::get_instance();
    //i must always be smaller than j
    //TODO abfangen wenn i == j
    if(i > j) {
        unsigned int k = j;
        j=i;
        i=k;
    }
    //same as in R mut -1 because index starts at 0
    unsigned int res = (c.L*(c.L-1)/2) - ((c.L-i+1)*((c.L-i+1)-1)/2) + j - i - 1;

    return res;
}

void FunctionalSequence::writeKdsToFile(const std::string& filename) {
    std::ofstream outfile(filename);

    if (outfile.good())
    {
            std::for_each(kds.cbegin(), kds.cend(), [&outfile](const auto& entry)
            {
                outfile << entry << '\n';
            });
    }
}
void FunctionalSequence::writeEpistasisToFile(const std::string& filename) {
    std::ofstream outfile(filename);

    if (outfile.good())
    {
        std::for_each(epistasis.cbegin(), epistasis.cend(), [&outfile](const auto& entry)
        {
            outfile << entry << '\n';
        });
    }
}

const double &FunctionalSequence::getEpistasis(const unsigned int i, const unsigned int j) {
    return epistasis[FunctionalSequence::getMatrixVectorIndex(i, j)];
}

const std::vector<double> &FunctionalSequence::getKd() {
    return kds;
}

const std::vector<double> &FunctionalSequence::getEpistasis() {
    return epistasis;
}
