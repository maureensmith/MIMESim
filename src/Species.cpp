//
//  Species.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 03.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Utils.hpp"
#include "Species.hpp"
#include "Constants.hpp"
#include "FunctionalSequence.hpp"
#include <numeric>
#include <algorithm>
#include <valarray>
#include <iostream>
#include <chrono>
#include <random>
#include <set>
#include <cmath>

namespace species
{
    Species::Species(const unsigned int id): specId(id), numMut(getNumberOfMutationsById()), mutatedPositions(specIdxToMutPos()), read(createRead()),
                                             count(0), mutCountBound(0), mutCountUnbound(0), errorCountBound(0.0), errorCountUnbound(0.0)
                                             {
                                             }

    unsigned int Species::getNumberOfMutationsById() {
        return species::getNumberOfMutationsById(this->specId);
    }

    mutVector Species::specIdxToMutPos() {
            return species::specIdxToMutPos(this->specId);
    }

    ref::ref_map Species::createRead() {
        constants::Constants& constants = constants::Constants::get_instance();
        ref::ref_map read;
        read.reserve(constants.L);
        for(unsigned i=1; i <= constants.L; ++i) {
            //TODO Workaround to set wildtype to A and mutations to C
            read.add({i, nucleotid::nucleobase{1}});
        }

        for(auto& mutPos:mutatedPositions) {
            auto pos = std::find_if(read.begin(), read.end(), [&mutPos](const auto& val)
            {
                return val.first == mutPos.getPosition();
            });
            if(pos != read.end())
                read.remove(pos);
            //TODO workaround: add the mutation as nucleobase, here simply  + 1 (->wt=A=1), create a wt to mut interpretation, or make it more general also for AA
            read.add({mutPos.getPosition(), nucleotid::nucleobase{int(mutPos.getSymbol()+1)}});
        }
        return read;
    }

    const unsigned int Species::getSpecId() const {
        return specId;
    }

    unsigned int Species::getCount() const {
        return count;
    }

    const unsigned int Species::getNumMut() const {
        return numMut;
    }

    const mutVector &Species::getMutatedPositions() const {
        return mutatedPositions;
    }

    double Species::getFreq() const {
        constants::Constants& c = constants::Constants::get_instance();
        return this->count/double(c.M);
    }

    double Species::getKd() const{
        return kd;
    }

    void Species::setCount(unsigned int count) {
        Species::count = count;
    }

    const ref::ref_map &Species::getRead() const {
        return read;
    }


    unsigned int Species::getMutCountBound() const {
        return mutCountBound;
    }

    void Species::setMutCountBound(unsigned int mutCountBound) {
        Species::mutCountBound = mutCountBound;
    }

    unsigned int Species::getMutCountUnbound() const {
        return mutCountUnbound;
    }

    void Species::setMutCountUnbound(unsigned int mutCountUnbound) {
        Species::mutCountUnbound = mutCountUnbound;
    }

    int Species::getErrorCountBound() const {
        return errorCountBound;
    }

    void Species::setErrorCountBound(int errorCountBound) {
        Species::errorCountBound = errorCountBound;
    }

    void Species::addErrorCountBound(int errorCountBound) {
        Species::errorCountBound += errorCountBound;
    }

    int Species::getErrorCountUnbound() const {
        return errorCountUnbound;
    }

    void Species::setErrorCountUnbound(int errorCountUnbound) {
        Species::errorCountUnbound = errorCountUnbound;
    }

    void Species::addErrorCountUnbound(int errorCountUnbound) {
        Species::errorCountUnbound += errorCountUnbound;
    }

    double Species::getTotalFractionBound(){
        constants::Constants& c = constants::Constants::get_instance();
        return mutCountBound/double(c.M);
    }

    double Species::getTotalFractionUnbound() {
        constants::Constants& c = constants::Constants::get_instance();
        return mutCountUnbound/double(c.M);
    }

    double Species::getFractionBound() {
        return getTotalFractionBound()/(getTotalFractionBound() + getTotalFractionUnbound());
    }

    double Species::getFractionUnbound() {
        return getTotalFractionUnbound()/(getTotalFractionBound() + getTotalFractionUnbound());
    }


    void Species::computeSpeciesKd() {
        FunctionalSequence& effects = FunctionalSequence::get_instance();
        constants::Constants& c = constants::Constants::get_instance();
        Species::kd = 1.0;
        //TODO nachfragen: wie setzt sich der Gesamteffekt zusammen? prod(Kd_i)*prod(e_ij) ?  oder prod(e_ij^2)?
        //additive effect of epistasis (since we have the exponential of the epistasis here, it is multiplicative
        for(auto mutPos1_it = begin(Species::mutatedPositions); mutPos1_it != end(Species::mutatedPositions) ; ++mutPos1_it) {
            Species::kd *= effects.getKd(*mutPos1_it);
            //std::cout << "mut pos1 " << *mutPos1_it << std::endl;
            for(auto mutPos2_it = mutPos1_it + 1; mutPos1_it != Species::mutatedPositions.end() && mutPos2_it != Species::mutatedPositions.end(); ++mutPos2_it) {
                //std::cout << "mut pos2 " << *mutPos2_it << std::endl;
                //compute the kd for a species by multiplying all single kds of the mutations and add the pairwise epistasis factor
                //totalEpistasisPerPos[*mutPos1_it-1] *= effects.getEpistasis(*mutPos1_it-1, *mutPos2_it-1);
                //totalEpistasisPerPos[*mutPos2_it-1] *= effects.getEpistasis(*mutPos1_it-1, *mutPos2_it-1);
                Species::kd *= effects.getEpistasis(*mutPos1_it, *mutPos2_it);
            }

        }
    }

    idCountMap drawSpeciesIds() {
        auto& constants = constants::Constants::get_instance();
        const auto seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator (seed);

        // Break down the drawing of all possible (allowed) species ids into 2 smaller ones:
        //first draw a the number of mutations from 0 to MAX_MUT, with the given probabilities...
        std::discrete_distribution<> d(begin(constants.P_NMUT), end(constants.P_NMUT));
        //then draw uniformly the id from the id range for this particular number of mutations
        std::vector<std::uniform_int_distribution<>> unif(constants.MAX_MUT);

        // create distributions for all numbers of mutations beforehand
        for(int numMut = 1; numMut<= constants.MAX_MUT; ++numMut) {
            unif[numMut-1] = std::uniform_int_distribution<>(constants.NMUT_RANGE[numMut - 1] + 1,
                                                             constants.NMUT_RANGE[numMut]);
        }
        // count the given species
        idCountMap m;
        for(int n=0; n<constants.M; ++n) {
            //draw number of mutations
            const int numMut = d(generator);
            //if no mutations, the id is always 1
            int id = 1;
            if(numMut > 0) {
                id = unif[numMut-1](generator);
                //std::cout << "id " << id << std::endl;
            }
            ++m[id];
        }
        return m;
    }

    //TODO löschen
//    idCountMap drawWildtypeErrors() {
//        auto& constants = constants::Constants::get_instance();
//        const int id = 1;
//        const int numMut = 0;
//        // count the given species
//        idCountMap m;
//        for(int n=0; n<constants.M; ++n) {
//            ++m[id][0];
//            //compute id with sequening errors:
//            auto newId = drawError(id, numMut);
//            --m[id][1];
//            ++m[newId][1];
//        }
//        return m;
//    }

    //TODO schneller machen in dem nicht immer neuer generator und verteilung erstellt wird?
    //TODO testen für q>2
    unsigned drawError(const unsigned id, const unsigned numMut) {
        auto& constants = constants::Constants::get_instance();
        const auto seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
        std::default_random_engine generator (seed);
        std::binomial_distribution<int> bino(constants.L,constants.P_ERR);
        //random generator for the position with an error
        std::uniform_int_distribution<> unif_err(1,constants.L);
        //random generator for the mutated symbol (if sequence symbol is wt, the errror is one of the mutations;
        // if the sequence symbol is mutated, the according symbol is the wild type or one of the other mutations
        std::uniform_int_distribution<> unif_sym(1,constants.Q-1);

        // if no errors: id stays the same
        unsigned newId = id;
        // draw number of errors, but only allow a maxmimum of MAX MUT mutations+errors in total //TODO vllt mal ändern und nicht nur bis 4 gehen?
        int numErrors = bino(generator);
        //TODO mal ausporbieren wenn ich mehr error erlaube,
        numErrors = std::min<int>(numErrors, constants.MAX_MUT*2 - numMut);

        if(numErrors > 0) {
            mutVector mutations = species::specIdxToMutPos(id);
            //containing numErrors positions with error
            //std::set<int> uniquePositions;
            //containing numErrors Mutations with unique position
            std::set<Mutation> uniquePositions;
            //TODO für q>2 anpasseninsert
            while(uniquePositions.size() < numErrors) {
                //draw position with error
                auto position = unif_err(generator);
                //...and the symbol
                auto symbol = unif_sym(generator);
                //first time we sample this position
                uniquePositions.emplace(position, symbol);
            }

            //if a real mutation has error, the according symbol need to be updated. In case it turns into wild type delete it
            for(auto& mut : mutations) {
                auto it=uniquePositions.find(mut);
                if(it != uniquePositions.end()) {
                    //TODO QUESTION teuer, aber selten.... machen?
                    //in case the real mutations symbol is drawn, read it as wild type and delete from list
                    if(it->getSymbol() == mut.getSymbol())
                        uniquePositions.erase(it);
                }
                else
                    uniquePositions.insert(mut);
            }

            //create a new species which includes both real and error mutations
            mutVector newMutations(uniquePositions.begin(), uniquePositions.end());
//            if(newMutations.size() > constants.MAX_MUT)
//                std::cout << "ACHTUNG: MEHR " << newMutations.size() << std::endl;

            newId = species::mutPosToSpecIdx(newMutations);
        }

        return newId;
    }

    //TODO: rekursiver Aufruf? Aber dafür müsste jedesmal für irgendein L' (Restlänger nach aktueller Position) die ID ranges berechnet werden, oder mache ich das eh?
    mutVector specIdxToMutPos(const unsigned specId) {
        constants::Constants& constants = constants::Constants::get_instance();
        auto numMut = getNumberOfMutationsById(specId);
        //collect the mutated position with the respective mutation symbol
        mutVector mutPos;
        mutPos.reserve(numMut);

        //TODO mit q anpassen
        //TODO weg: save the "symbol" of the mutation for each position
        //mutVector mutSymbol(numMut);
        // check if Id is valid
        if(specId<=constants.NMUT_RANGE.back() && numMut > 0) {
            unsigned int Lact = constants.L;
            unsigned int numMutAct = numMut;
            auto mSymbols = constants.Q-1;

            //TODO weg
//            if(specId ==22)
//                std::cout << "Blub";
            //get the id within the range of number of mutations (substract the ids for the sequences with less mutations)
            unsigned int idAct = specId - constants.NMUT_RANGE[numMut-1];


            // determine each mutations position seen from the mutations position before...
            for(unsigned int m=0; m<numMut; ++m){
                // for each possible positions within the length the actual mutation covers a range of ids depending on the residual mutations to follow
                std::vector<unsigned long> cumSumRange(Lact-(numMutAct-1));
                //bool indexFound = false;
                unsigned int i = 0;
                //initialise first value of the vector for cummulative sum (do it so complicated to not compute the whole range if not necessary)
                cumSumRange[0] = utils::nChoosek(Lact-1, numMutAct-1) * std::pow(mSymbols, numMut);
                // find the id within the ranges and get the index (=position)
                while(idAct > cumSumRange[i] && i<Lact-(numMutAct-1)) {
                    ++i;
                //for(; i<Lact-(numMutAct-1) && !indexFound; ++i) {
                    //TODO test
                    cumSumRange[i] = utils::nChoosek(Lact-i-1, numMutAct-1) * std::pow(mSymbols, numMut);
                    if(i>0) {
                        cumSumRange[i] += cumSumRange[i-1];
                    }
                }

                //TODO weg
                //mutPos[m] = i + 1;
                // the symbol of a posisition  is given for (q-1) ^ numMut-1 times (e.g. with three mutations and 2 symbols 2 ^2 times) : AAA, AAB, ABA, ABB
                int symbolCombiPerPos = std::pow(mSymbols, numMutAct-1);
                //TODO getestet? wenn ja, comments weg
                //find symbol
                unsigned mut = (int)std::floor((idAct-1) / symbolCombiPerPos) % mSymbols;

                unsigned int pos = i + 1; //TODO Utils hat noch das- *cumSumRange.begin()
                unsigned int prePos = mutPos.begin() == mutPos.end() ? 0 : mutPos.rbegin()->getPosition();

                //arguments: the two pair_constructor parameter pos and mut, adding the last cummulative position
                // (= position seend frim the beginning of the sequence)
                mutPos.emplace_back(pos + prePos, mut);

                       // indexFound = true;

                //the residual length after the actual mutation
                Lact = Lact - pos;
                // the redsiudal number of mutations after the actual mutation
                --numMutAct;
                // the id within the residual length (-1 because the indices start at 0 and another -1 because we
                // substract the ids of the preceeding mutations range
                idAct =  idAct - (pos==1 ? 0 : cumSumRange[i-1]);// hier stimmt was nicht

            }
        }
        //.... and get the correct positions within the sequence with cumsum
        //TODO in die schleife rein, wie in Utils..
        //std::partial_sum(mutPos.begin(), mutPos.end(), mutPos.begin());

        //TODO weg wenn alles getestet
        //if(numMut==5) {
         //   std::cout << "id: " << specId << std::endl;
//            for (auto blub : mutPos) {
//                std::cout << "pos " << blub.getPosition() << " mut " << blub.getSymbol() << std::endl;
//            }
    //    }
            //TODO weg, hab ich ja shcon oben gelöst
//            std::partial_sum(mutPos.begin(), mutPos.end(), mutPos.begin(),
//                             [](const Mutation &x, const Mutation &y) {
//                                 return Mutation(x.getPosition() + y.getPosition(), y.getSymbol());
//                             }
//            );
//            std::cout << "nachher: " << std::endl;
//            for (auto blub : mutPos) {
//                std::cout << "pos " << blub.getPosition() << " mut " << blub.getSymbol() << std::endl;
//            }
//        }

        return(mutPos);
    }

    //TODO testen! (vorallem das mit partial sum)
    unsigned mutPosToSpecIdx(const mutVector& mutPos) {
        constants::Constants& constants = constants::Constants::get_instance();
        unsigned numMut = mutPos.size();
        //id for 0 mutations is 1
        unsigned specId = 1;
        if(numMut > 0) {
            // add the ids for the sequences with less mutations
            specId = constants.NMUT_RANGE[numMut-1];
            mutVector mutPos_new = mutPos;
            // get the indices for the individual length segments for each position
            std::partial_sum(mutPos.begin(), mutPos.end(), mutPos_new.begin(),
                    [](const Mutation & x, const Mutation& y){return Mutation(y.getPosition() - x.getPosition(), y.getSymbol());});
            //TODO weg get the indices for the individual length segments for each position
//             for(unsigned i = 1; i<mutPos_new.size(); ++i) {
//                 mutPos_new[i].getPosition() -= mutPos[i-1].getPosition();
//             }
            unsigned Lact = constants.L;
            unsigned numMutAct = numMut;

            //Notiz an mich selbst: enforcing const elements in range iteration (C++17)
            for(auto const& mutation : std::as_const(mutPos_new)) {
                if(numMutAct == 1) {
                    //TODO überall checken: the possible mutations are represented from 0 upwards
                    specId += mutation.getPosition()+mutation.getSymbol();
                } else {
                    // if the position of the actual mutations is != 1, the id is depending on the next position
                    if( mutation.getPosition()!=1) {
                        for(int i = 1; i<mutation.getPosition(); ++i) {
                            // specId += utils::nChoosek(Lact-i, numMutAct-1);
                            specId += utils::nChoosek(Lact-i, numMutAct-1) * pow(constants.Q-1, numMutAct);
                        }
                    }
                    Lact -= mutation.getPosition();
                    --numMutAct;
                }
            }
        }

        return specId;
    }

    unsigned getNumberOfMutationsById(const unsigned specId) {
        constants::Constants& constants = constants::Constants::get_instance();
        //gives the index where the content is still lower than the given id
        auto low_it = std::lower_bound(std::begin(constants.NMUT_RANGE), std::end(constants.NMUT_RANGE), specId);
        return (low_it - std::begin(constants.NMUT_RANGE));
    }

    //TODO: Umbau nach counts
    //void addCountsWithError(const std::valarray<double>& freqBound, const std::valarray<double>& freqUnbound, species_map& spec_map) {
    void addCountsWithError(const std::valarray<unsigned int>& SBound, const std::valarray<unsigned int>& SUnbound, species_map& spec_map) {
        constants::Constants& cons = constants::Constants::get_instance();
        // collecting the bound(0) and unbound(1) counts for species which occur due to noise
        std::map<int, std::array<int,2>> errCounts;
        // to get the correct counts from the valarrays
        int specIdx = 0;
        for(auto it = spec_map.begin(); it != spec_map.end(); ++it) {
            unsigned oldId = it->first;
            // set the counts
            //TODO: Umbau nach counts
            //it->second.setMutCountUnbound(std::round(freqUnbound[specIdx]*cons.M));
            //it->second.setMutCountBound(std::round(freqBound[specIdx]*cons.M));
            it->second.setMutCountUnbound(SUnbound[specIdx]);
            it->second.setMutCountBound(SBound[specIdx]);

            // sample error for bound sequences
            for(int b = 0; b<(it->second.getMutCountBound());++b) {
                unsigned long newId = drawError(oldId, it->second.getNumMut());
                auto newIt = spec_map.find(newId);
                if (newIt != spec_map.end()) {
                    if(oldId != newId) {
                        newIt->second.addErrorCountBound(1);
                        it->second.addErrorCountBound(-1);
                    }
                } else {
                    // in case a species is not sampled yet, put into the map with its error counts later.
                    ++errCounts[newId][0];
                    it->second.addErrorCountBound(-1);
                }
            }
            //sample error for unbound sequences
            for(int b = 0; b<(it->second.getMutCountUnbound());++b) {
                unsigned long newId = drawError(oldId, it->second.getNumMut());
                auto newIt = spec_map.find(newId);
                if (newIt  != spec_map.end()) {
                    if(oldId != newId) {
                        newIt->second.addErrorCountUnbound(1);
                        it->second.addErrorCountUnbound(-1);
                    }

                } else {
                    // in case a species is not sampled yet, put into the map with its error counts later.
                    ++errCounts[newId][1];
                    it->second.addErrorCountUnbound(-1);
                }
            }
            ++specIdx;
        }

        // add species which only occur because of sequencing errors
        for(auto it = errCounts.begin(); it != errCounts.end(); ++it) {
            auto currentObj = spec_map.emplace(it->first, it->first);
            // first is a pointer to just constructed pair, second is the species object to call the methods
            currentObj.first->second.setErrorCountBound(it->second[0]);
            currentObj.first->second.setErrorCountUnbound(it->second[1]);
        }

    }

}
