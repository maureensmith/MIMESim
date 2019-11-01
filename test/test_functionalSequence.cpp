//
// Created by Smith, Maureen on 31.05.18.
//


#include <iostream>
#include "catch.hpp"
#include "FunctionalSequence.hpp"
#include "Constants.hpp"

TEST_CASE("Testing FunctionalSequence Class") {
    //TODO: wie kann ich für verschiedene Test cases verschieden set ups erstellen so dass sie innerhalb des cases für alle sections gilt, aber nicht für alle cases?
    const unsigned int length = 10;
    const unsigned int q = 2;
    const double p_mut = 0.1;
    const double p_error = p_mut/10;
    const double p_effect = 0.5;
    const double p_epistasis = 0.3;
    constants::Constants &cons = constants::Constants::create_instance(length, q, p_mut, p_error, p_effect,  p_epistasis);

    SECTION("test function getMatrixIndex which converts 2 indicices into the corresponding index  of a vectoral symmetric matrix") {

        FunctionalSequence& effect = FunctionalSequence::get_instance();
        unsigned int i;
        unsigned int j;
        unsigned int idx = 0;

        auto epistasis = effect.getEpistasis();

        for(i=1; i<cons.L; ++i) {
            for(j=i+1; j<=cons.L; ++j) {
                REQUIRE(effect.getEpistasis(i, j) == epistasis[idx]);
                ++idx;
            }
        }
    }

    SECTION("Test if its singleton") {
        FunctionalSequence& effect = FunctionalSequence::get_instance();
        // Kds are picked randomly. If it would not be a singleton, the second drawn effects would be different
        FunctionalSequence& effect2 = FunctionalSequence::get_instance();
        for(int i = 1; i<=cons.L; ++i) {
            REQUIRE(effect.getKd(i) == effect2.getKd(i));
        }
    }

    SECTION("Test total species Kd computation") {
        //TODO Test schreiben

    }


}