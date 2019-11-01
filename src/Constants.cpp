//
//  Constants.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 04.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Constants.hpp"
#include <exception>
#include <math.h>
#include <iostream>

namespace constants {
    // separate definitions for static constexpr need to be placed out of the class to enable linkage
    constexpr double Constants::KD_WT;
    //constexpr double Constants::P_EFFECT;
    constexpr double Constants::NO_EPISTASIS;
    //constexpr double Constants::P_EPISTASIS;
    //constexpr unsigned int Constants::MAX_MUT ;
    constexpr unsigned int Constants::M;

    /*
    * determination of the maximal number of mutations (Expected number of occurrence >5), must be minimum 3
    */
    unsigned int Constants::computeMaxMut(const unsigned int L, const double pMut) {
        unsigned int maxMut = 2;
        do{
            ++maxMut;
            std::cout << "pmut " << pMut << " L " << L << std::endl;
            std::cout << "max mut" << maxMut << " Exp. count " << M * pow(pMut,maxMut+1) * pow(1-pMut, L-maxMut-1) << std::endl;
        }while(M * pow(pMut,maxMut+1) * pow(1-pMut, L-maxMut-1) > 5 && maxMut < L);
        return maxMut;
    }

    //std::array<unsigned int, Constants::MAX_MUT+1> Constants::setNMutRange() {
    std::vector<unsigned int> Constants::setNMutRange(const unsigned int maxMut, const unsigned int L) {
        //compute the number of possible sequence for 0..MAX_MUT mutations, the cumulative sum gives the id range for each number of mutations
        //TODO mal ausporbieren wenn ich mehr error erlaube, also einfach nur die größere range für spätere eror id berechnung
        std::vector<unsigned int> nMutRange(maxMut*2+1);
        //std::array<unsigned int, MAX_MUT+1> nMutRange;
        nMutRange[0] = 1;
        for(unsigned int i = 1; i<=maxMut*2; ++i) {
            nMutRange[i] = nMutRange[i - 1] + utils::nChoosek(L, i);
        }
        return(nMutRange);
    }

    /*
     * computing the probability for each number of mutations
     */
    //std::array<double, Constants::MAX_MUT+1> Constants::setP_NMut() {
    std::vector<double> Constants::setP_NMut(const unsigned int maxMut, const unsigned int L, const double pMut) {
        //std::array<double, Constants::MAX_MUT+1> p_nmut;
        std::vector<double> p_nmut(maxMut+1);
        double p_sum = 0;
        //TODO test: die wahrscheinlichekiten für n > n_max werden zu n_max hinzugezählt
        //for(unsigned i = 0; i<= maxMut; ++i) {
        for(unsigned i = 0; i< maxMut; ++i) {
            p_nmut[i] = utils::nChoosek(L, i)* pow(pMut,i)*pow(1-pMut,L-i);
            p_sum += p_nmut[i];
        }
        p_nmut[maxMut] = 1-p_sum;
        return(p_nmut);
    };

    Constants &Constants::get_instance() {
        constexpr unsigned int L_mock = 0;
        Constants& c = create_instance(L_mock, 0, 0, 0, 0, 0);
        if(L_mock != c.L)
            return c;
        else
            throw("Singleton not correctly instantiated.");
    }

//    Constants& Constants::create_instance(unsigned int length, unsigned int q, double p_mut){
//        return create_instance(length, q, p_mut, 0);
//    }

    Constants& Constants::create_instance(unsigned int length, unsigned int q, double p_mut, double p_error, double p_effect, double p_epistasis){
        //static Constants instance(length);
        static Constants instance(length, q, p_mut, p_error, p_effect, p_epistasis);
        return instance;
    }
}