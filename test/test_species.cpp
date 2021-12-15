//
//  test_species.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 24.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Constants.hpp"
#include "Species.hpp"
#include "catch.hpp"

#include <iostream>
#include <valarray>

void testSpecIdxToMutPos()
{
    constants::Constants& cons = constants::Constants::get_instance();
    auto q = cons.Q;
    auto length = cons.L;
    auto numSymbols = q - 1;

    // Test id 1 = no  mutation
    unsigned int id = 1;
    auto positions = species::specIdxToMutPos(id);
    REQUIRE(positions.size() == 0);
    REQUIRE(positions.empty());

    // test Id range for one mutation
    for (auto pos = 1; pos <= length; ++pos)
    {
        for (auto symb = 1; symb < q; ++symb)
        {
            ++id;
            positions = species::specIdxToMutPos(id);
            REQUIRE(positions.size() == 1);
            REQUIRE(positions.back() == pos);
        }
    }

    for (auto pos1 = 1; pos1 <= length - 1; ++pos1)
    {
        for (auto pos2 = pos1 + 1; pos2 <= length; ++pos2)
        {
            for (auto symb1 = 1; symb1 < q; ++symb1)
            {
                // for (auto symb2 = 1; symb2 < q; ++symb2) {
                ++id;
                positions = species::specIdxToMutPos(id);
                REQUIRE(positions.size() == 2);
                REQUIRE(positions.at(0) == pos1);
                REQUIRE(positions.at(1) == pos2);
                //}
            }
        }
    }

    //    for (auto pos1 = 1; pos1 <= length - 2; ++pos1) {
    //        for (auto pos2 = pos1 + 1; pos2 <= length-1; ++pos2) {
    //            for (auto pos3 = pos2 + 1; pos3 <= length; ++pos3) {
    //            for (auto symb = 1; symb < q; ++symb) {
    //                ++id;
    //                positions = species::specIdxToMutPos(id);
    //                REQUIRE(positions.size() == 2);
    //                REQUIRE(positions.at(0) == pos1);
    //                REQUIRE(positions.at(1) == pos2);
    //            }
    //        }
    //
    //    }
}

TEST_CASE("Testing Species Utility methods")
{

    constants::Constants& cons = constants::Constants::get_instance();

    SECTION("Test function drawSpeciesId")
    {
        species::idCountMap idCounts = species::drawSpeciesIds();
        std::valarray<int> numMutCounts(cons.MAX_MUT + 1);
        for (auto it = idCounts.begin(); it != idCounts.end(); ++it)
        {
            // Test if the id range is valid
            REQUIRE(it->first > 0);
            REQUIRE(it->first <= cons.NMUT_RANGE[cons.MAX_MUT]);
            species::Species spec(it->first);
            numMutCounts[spec.getNumMut()] += it->second;
        }
        // Test that all number of mutations were drawn
        for (auto& count : numMutCounts)
        {
            REQUIRE(count > 0);
        }
        REQUIRE(constants::Constants::M == numMutCounts.sum());
    }

    SECTION("Test function specIdxToMutPos")
    {
        testSpecIdxToMutPos();
    }

    SECTION("Test function mutPosToSpecIdx")
    {
        for (int i = 1; i <= cons.NMUT_RANGE[cons.MAX_MUT]; ++i)
        {
            unsigned id = i;
            auto mutPos = species::specIdxToMutPos(id);
            unsigned computedId = species::mutPosToSpecIdx(mutPos);
            REQUIRE(id == computedId);
        }
    }

    SECTION("Test function addCountsWithError")
    {
        std::valarray<double> f_bound_tot;
        std::valarray<double> f_unbound_tot;
    }
}

TEST_CASE("Testing Species Utility methods with q > 2")
{
    // TODO: wie kann ich für verschiedene Test cases verschieden set ups erstellen so dass sie innerhalb des cases für
    // alle sections gilt, aber nicht für alle cases?
    const unsigned int length = 10;
    const unsigned int q = 4;
    const double p_mut = 0.1;
    const double p_error = p_mut / 10;
    const double p_effect = 0.5;
    const double p_epistasis = 0.3;
    constants::Constants& cons = constants::Constants::create_instance(length, q, p_mut, p_error, p_effect, p_epistasis,
                                                                       std::filesystem::path());
    std::cout << "Q " << cons.Q << std::endl;
    SECTION("Test function specIdxToMutPos")
    {
        testSpecIdxToMutPos();
    }
}

TEST_CASE("Testing Species Class")
{

    constants::Constants& cons = constants::Constants::get_instance();
    // 12*10^6 * 0.1^6 * 0.9^4 = 7.8732
    REQUIRE(cons.MAX_MUT == 6);
    // TODO anpassen, wenn entschieden ist, wieviel erros (=mutationen) erlaubt sind
    // REQUIRE(cons.NMUT_RANGE.size() == cons.MAX_MUT*2+1);

    SECTION("test function getNumberOfMutationsById")
    {

        unsigned int numMut = 2;
        unsigned int id = 12;
        species::Species spec(id);
        REQUIRE(spec.getNumMut() == numMut);

        numMut = 0;
        id = 1;
        species::Species spec2(id);
        REQUIRE(spec2.getNumMut() == numMut);
    }

    SECTION("test function idToMutatedPositions")
    {
        unsigned int numMut = 2;
        unsigned int id = 12;
        species::Species spec(id);
        species::mutVector mutPos = spec.getMutatedPositions();
        // std::cout << "length " << mutPos.size() << std::endl;
        REQUIRE(mutPos.size() == numMut);
        REQUIRE(mutPos[0] == 1);
        REQUIRE(mutPos[1] == 2);
    }

    SECTION("test function createRead")
    {
        unsigned int id = 12;
        species::Species spec(id);
        auto& read = spec.getRead();
        // TODO test if the read has the mutations at the correct position and length is L

        unsigned size = read.size();
        REQUIRE(size == cons.L);
    }
}
