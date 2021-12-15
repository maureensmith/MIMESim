//
//  Species.hpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 03.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#ifndef Species_hpp
#define Species_hpp

#include "Count.hpp"
#include "Mutation.hpp"

#include <array>
#include <map>
#include <random>
#include <set>
#include <stdio.h>
#include <unordered_map>
#include <vector>

namespace species
{
    // TODO weg nach test
    // typedef std::array<unsigned int, Constants::MAX_MUT> mutArr;
    // typedef std::vector<unsigned int> mutVector;
    // TODO ist ja eigentlich generell, nicht nur Species. Irgedwo anders hinschieben?
    // TODO stattdessen einfach Mutation?
    // typedef std::pair<unsigned, unsigned> posMutPair;
    // contain each mutated position with the respective mutation
    typedef std::vector<Mutation> mutVector;
    // typedef std::map<int, std::array<int, 2>> idCountMap;
    // TODO umändern in unordered_map und vergleichen
    typedef std::unordered_map<int, int> idCountMap;

    class Species
    {

      private:
        const unsigned int specId;
        unsigned int count;
        const unsigned int numMut;
        // mutatated positions need to be in ascending order
        // TODO where to test this, exception?
        const mutVector mutatedPositions;
        // double freq;
        // KD of the the given sequence, adding all single Kds of the mutations and the epistatic effects for pairs
        double kd = 1.0;
        // TODO for later count purposes directly create a read. nee weg damit, nimmt zu viel platz ein
        // const ref::ref_map read;

        /**** after ODE: save bound and unbound fraction information ****/

        unsigned mutCountBound;
        unsigned mutCountUnbound;

        int errorCountBound;
        int errorCountUnbound;

        unsigned int getNumberOfMutationsById();

        mutVector specIdxToMutPos();

      public:
        Species(const unsigned int id);

        const unsigned int getSpecId() const;

        unsigned int getCount() const;

        const unsigned int getNumMut() const;

        const mutVector& getMutatedPositions() const;

        double getFreq() const;

        double getKd() const;

        void setCount(unsigned int count);

        void incrementCount();

        void computeSpeciesKd();

        // TODO weg
        // const ref::ref_map &getRead() const;
        // ref::ref_map createRead();

        unsigned int getMutCountBound() const;

        void setMutCountBound(unsigned int mutCountBound);

        unsigned int getMutCountUnbound() const;

        void setMutCountUnbound(unsigned int mutCountUnbound);

        int getErrorCountBound() const;

        void setErrorCountBound(int errorCountBound);

        void addErrorCountBound(int errorCountBound);

        int getErrorCountUnbound() const;

        void setErrorCountUnbound(int errorCountUnbound);

        void addErrorCountUnbound(int errorCountUnbound);

        // Fraction of all M sequences
        double getTotalFractionBound();

        double getTotalFractionUnbound();

        // fraction of the one particular species
        double getFractionBound();

        double getFractionUnbound();
    };

    using species_map = std::unordered_map<int, Species>;

    mutVector specIdxToMutPos(const unsigned specId);

    unsigned mutPosToSpecIdx(const mutVector& mutPos);

    unsigned getNumberOfMutationsById(const unsigned specId);

    species_map drawSpeciesIds();

    std::set<Mutation> drawError_2(const mutVector& mutations, std::default_random_engine& generator);

    count::counter_collection countMutationsWithErrors(const std::valarray<unsigned int>& SBound,
                                                       const std::valarray<unsigned int>& SUnbound,
                                                       const species_map& spec_map);

}
#endif /* Species_hpp */
