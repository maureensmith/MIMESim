//
// Created by Smith, Maureen on 07.06.18.
//

#include <iostream>
#include <valarray>
#include "catch.hpp"
#include "Utils.hpp"
#include "Constants.hpp"
#include "Species.hpp"

TEST_CASE("Testing Utils") {
    //TODO: wie kann ich für verschiedene Test cases verschieden set ups erstellen so dass sie innerhalb des cases für alle sections gilt, aber nicht für alle cases?
    const unsigned int length = 10;
    const unsigned int q = 2;
    const double p_mut = 0.1;
    const double p_error = p_mut/10;
    const double p_effect = 0.5;
    const double p_epistasis = 0.3;
    constants::Constants &cons = constants::Constants::create_instance(length, q, p_mut, p_error, p_effect,  p_epistasis);

    SECTION("Test  function n choose k (binomial coefficient)") {

    }

}