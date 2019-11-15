//
//  Constants.hpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 04.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#ifndef Constants_hpp
#define Constants_hpp

#include <stdio.h>
#include <array>
#include "Utils.hpp"
#include <math.h>

namespace constants {
    struct Constants {
        //static = lifetime during execution (run time);
        //constexpr value available at compile time -> constant,  (not necessarily run time)

        /***** Constants regarding output *******/
        //directory where the results and parameter file are stored
        static const std::string OUTPUT_DIR;
        //parameter file name never changes
        static constexpr auto PARAMETER_FILE = "parameters.txt";

        /***** Constants regarding kds sampling *******/
        //absolute wildtype Kd
        static constexpr double KD_WT = 1.0;
        //probability for each position to have an effect when mutated
        //static constexpr double P_EFFECT = 0.5;
        const double P_EFFECT = 0.5;
        //no multiplicative effect, when there is no epistasis
        static constexpr double NO_EPISTASIS = 1.0;
        //probability for each position pair to have epistatic effects
        //static constexpr double P_EPISTASIS = 0.75;
        const double P_EPISTASIS = 0.75;


        /**** Constants regarding species sampling *******/
        // maximal number of mutations in one sequence
        //static constexpr unsigned int MAX_MUT = 3;
        const unsigned int MAX_MUT;
        //number of total sequences
        static constexpr unsigned int M = 12*pow(10,6);
        //sequence length
        const unsigned int L;
        //number of pairwise values
        const unsigned int PW;
        // symbols per position
        const unsigned int Q;
        // probability for a mutation
        const double P_MUT = 0.01;
        //probability for sequencing error
        const double P_ERR = 0.001;
        // id range for 0..MAX_MUT mutations for a sequence length L: (0 -> 1, 1 -> 2..L+1, etc)
        //const std::array<unsigned int, MAX_MUT+1> NMUT_RANGE;
        const std::vector<unsigned int> NMUT_RANGE;
        //const std::array<double, MAX_MUT+1> P_NMUT;
        const std::vector<double> P_NMUT;

        //std::array<unsigned int, MAX_MUT+1> setNMutRange();
        //std::array<double, MAX_MUT+1> setP_NMut();
        std::vector<unsigned int> setNMutRange(const unsigned int maxMut, const unsigned L );
        std::vector<double> setP_NMut(const unsigned int MaxMut, const unsigned L, const double pMut);
        unsigned int computeMaxMut(const unsigned int L, const double pMut);

        static Constants& create_instance(const unsigned int length, const unsigned int q, const double p_mut,
                                          const double p_error, const double p_effect, const double p_epistasis);
        //static Constants& create_instance(const unsigned int length, const unsigned int q, const double p_mut);

        static Constants& get_instance();

        // The copy constructor is deleted, to prevent client code from creating new
        // instances of this class by copying the instance returned by create_instance()
        Constants(Constants const&) = delete;

        // The move constructor is deleted, to prevent client code from moving from
        // the object returned by create_instance(), which could result in other clients
        // retrieving a reference to an object with unspecified state.
        Constants(Constants&&) = delete;

    private:
        // Default-constructor is private, to prevent client code from creating new
        // instances of this class. The only instance shall be retrieved through the
        // create_instance() function.
        //TODO: workaround, 3 mal die gleiche routine (computeMaxMut) aufrufen, weil es erst am  ende alles gespeichert wird, anders lösen?
        Constants(unsigned int length, unsigned int q, double p_mut, double p_error, double p_effect, double p_epistasis) :
                L(length),PW(length*(length-1)/2), Q(q), P_MUT(p_mut), P_ERR(p_error), P_EFFECT(p_effect), P_EPISTASIS(p_epistasis),
                MAX_MUT(computeMaxMut(length,p_mut)), NMUT_RANGE(setNMutRange(computeMaxMut(length,p_mut),L)), P_NMUT(setP_NMut(computeMaxMut(length,p_mut),length,p_mut)) {};
        //Constants(unsigned int length, unsigned int q, double p_mut) : L(length),PW(L*(L-1)/2), Q(q), P_MUT(p_mut), NMUT_RANGE(setNMutRange()), P_NMUT(setP_NMut()) {};
    };
}
#endif /* Constants_hpp */
