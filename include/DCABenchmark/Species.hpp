//
//  Species.hpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 03.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#ifndef Species_hpp
#define Species_hpp

#include <stdio.h>
#include <array>
#include <vector>
#include <map>
#include "Mutation.hpp"

#include "sam2counts/ref_map.hpp"

namespace species {
    //typedef std::array<unsigned int, Constants::MAX_MUT> mutArr;
    typedef std::vector<unsigned int> posVector;
    //typedef std::map<int, std::array<int, 2>> idCountMap;
    typedef std::map<int, int> idCountMap;

    idCountMap drawSpeciesIds();

    unsigned drawError(const unsigned id, const unsigned numMut);

    posVector specIdxToMutPos(const unsigned specId);

    unsigned mutPosToSpecIdx(const posVector mutPos);

    unsigned getNumberOfMutationsById(const unsigned specId);


    
    class Species {

    private:
        const unsigned int specId;
        unsigned int count;
        const unsigned int numMut;
        //const unsigned int
        const posVector mutatedPositions;
        //double freq;
        double kd = 1.0;
        //for later count purposes directly create a read
        const ref::ref_map read;

        /**** after ODE: save bound and unbound fraction information ****/

        unsigned mutCountBound;
        unsigned mutCountUnbound;

        int errorCountBound;
        int errorCountUnbound;


        unsigned int getNumberOfMutationsById();

        posVector specIdxToMutPos();

        ref::ref_map createRead();
    public:
        Species(const unsigned int id);

        const unsigned int getSpecId() const;

        unsigned int getCount() const;

        const unsigned int getNumMut() const;

        const posVector &getMutatedPositions() const;

        double getFreq() const;

        double getKd() const;

        void setCount(unsigned int count);

        void computeSpeciesKd();

        const ref::ref_map &getRead() const;

        unsigned int getMutCountBound() const;

        void setMutCountBound(unsigned int mutCountBound);

        unsigned int getMutCountUnbound() const;

        void setMutCountUnbound(unsigned int mutCountUnbound);

        int getErrorCountBound() const;

        void setErrorCountBound(int errorCuntBound);

        void addErrorCountBound(int errorCuntBound);

        int getErrorCountUnbound() const;

        void setErrorCountUnbound(int errorCountUnbound);

        void addErrorCountUnbound(int errorCountUnbound);

        //Fraction of all M sequences
        double getTotalFractionBound();

        double getTotalFractionUnbound();

        //fraction of the one particular species
        double getFractionBound();

        double getFractionUnbound();

    };

    using species_map = std::map<int, Species>;

    //TODO: Umbau nach counts
    //void addCountsWithError(const std::valarray<double>& freqBound, const std::valarray<double>& freqUnbound,
     //                       species_map& spec_map);
    void addCountsWithError(const std::valarray<unsigned int>& SBound, const std::valarray<unsigned int>& SUnbound,
                            species_map& spec_map);

}
#endif /* Species_hpp */
