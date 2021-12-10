//
//  Constants.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 04.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Constants.hpp"

#include <exception>
#include <iostream>
#include <math.h>

namespace constants
{

    // Initialize static member inited
    bool Constants::inited = false;

    // separate definitions for static constexpr need to be placed out of the class to enable linkage
    constexpr double Constants::KD_WT;
    // constexpr double Constants::P_EFFECT;
    constexpr double Constants::NO_EPISTASIS;
    // constexpr double Constants::P_EPISTASIS;
    // constexpr unsigned int Constants::MAX_MUT ;
    constexpr unsigned int Constants::M;

    /*
     * determination of the maximal number of mutations (Expected number of occurrence >5), must be minimum 3
     */
    unsigned int Constants::computeMaxMut(const unsigned int L, const double pMut)
    {
        unsigned int maxMut = 2;
        do
        {
            ++maxMut;
        } while (M * pow(pMut, maxMut + 1) * pow(1 - pMut, L - maxMut - 1) > 5 && maxMut < L);
        return maxMut;
    }

    // TODO löschen
    //    //std::array<unsigned int, Constants::MAX_MUT+1> Constants::setNMutRange() {
    //    std::vector<unsigned int> Constants::setNMutRange(const unsigned int maxMut, const unsigned int L) {
    //        //compute the number of possible sequence for 0..MAX_MUT mutations, the cumulative sum gives the id range
    //        for each number of mutations
    //        //TODO mal ausporbieren wenn ich mehr error erlaube, also einfach nur die größere range für spätere eror
    //        id berechnung std::vector<unsigned int> nMutRange(maxMut*2+1);
    //        //std::array<unsigned int, MAX_MUT+1> nMutRange;
    //        nMutRange[0] = 1;
    //        for(unsigned int i = 1; i<=maxMut*2; ++i) {
    //            nMutRange[i] = nMutRange[i - 1] + utils::nChoosek(L, i);
    //        }
    //        return(nMutRange);
    //    }

    std::vector<unsigned int> Constants::setNMutRange(const unsigned int maxMut, const unsigned int L,
                                                      const unsigned int q)
    {
        // compute the number of possible sequence for 0..MAX_MUT mutations, the cumulative sum gives the id range for
        // each number of mutations larger range for adding errors (i.e. more mutated positions)
        std::vector<unsigned int> nMutRange(maxMut * 2 + 1);
        nMutRange[0] = 1;
        for (unsigned int i = 1; i <= maxMut * 2; ++i)
        {
            // cummulative sum of L choose i position combinations and (q-1)^i symbol combinations
            nMutRange[i] = nMutRange[i - 1] + utils::nChoosek(L, i) * pow(q - 1, i);
        }
        return (nMutRange);
    }

    /*
     * computing the probability for each number of mutations
     */
    // std::array<double, Constants::MAX_MUT+1> Constants::setP_NMut() {
    std::vector<double> Constants::setP_NMut(const unsigned int maxMut, const unsigned int L, const double pMut)
    {
        // std::array<double, Constants::MAX_MUT+1> p_nmut;
        std::vector<double> p_nmut(maxMut + 1);
        double p_sum = 0;
        // TODO test: die wahrscheinlichekiten für n > n_max werden zu n_max hinzugezählt
        // for(unsigned i = 0; i<= maxMut; ++i) {
        for (unsigned i = 0; i < maxMut; ++i)
        {
            p_nmut[i] = utils::nChoosek(L, i) * pow(pMut, i) * pow(1 - pMut, L - i);
            p_sum += p_nmut[i];
        }
        // the probabilities for n > n_max are added to n_max
        p_nmut[maxMut] = 1 - p_sum;
        return (p_nmut);
    };

    Constants &Constants::get_instance()
    {
        if (inited)
        {
            // Guaranteed to be destroyed.
            // Instantiated on first use.
            // mock call of constructor
            static Constants &c = create_instance(0, 0, 0, 0, 0, 0, std::filesystem::path());
            return c;
        }
        // constexpr unsigned int L_mock = 0;
        // Constants& c = create_instance(L_mock, 0, 0, 0, 0, 0, std::filesystem::path());
        // if(L_mock != c.L)
        //     return c;
        else
        {
            std::cerr << "Singleton not correctly instantiated." << std::endl;
            // TODO anderen exception type
            throw std::invalid_argument("Singleton not correctly instantiated.");
        }
    }

    Constants &Constants::create_instance(const unsigned int length, const unsigned int q, const double p_mut,
                                          const double p_error, const double p_effect, const double p_epistasis,
                                          const fs::path outputDir)
    {
        // static Constants instance(length);
        static Constants instance(length, q, p_mut, p_error, p_effect, p_epistasis, outputDir);
        inited = true;
        return instance;
    }
}
