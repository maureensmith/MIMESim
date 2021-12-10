//
// Created by msmith on 04.12.19.
//
#include "Constants.hpp"
#include "catch.hpp"

#include <iostream>

TEST_CASE("Testing Constants Class", "[ConstantClassTest]")
{
    SECTION("test get_instance, if there is no singleton instance yet")
    {
        // if no instance has been build before, an exception should be thrown.
        REQUIRE_THROWS_AS(constants::Constants::get_instance(), std::exception);
    }

    SECTION("test creation of constant class")
    {

        const unsigned int length = 10;
        const unsigned int q = 2;
        const double p_mut = 0.1;
        const double p_error = p_mut / 10;
        const double p_effect = 0.5;
        const double p_epistasis = 0.3;
        constants::Constants &cons = constants::Constants::create_instance(length, q, p_mut, p_error, p_effect,
                                                                           p_epistasis, std::filesystem::path());
        REQUIRE(cons.L == length);
        REQUIRE(cons.Q == q);
        REQUIRE(cons.P_MUT == p_mut);
        REQUIRE(cons.P_EFFECT == p_effect);
        REQUIRE(cons.P_EPISTASIS == p_epistasis);
        REQUIRE(cons.OUTPUT_DIR.empty());
    }

    SECTION("test creation of constant class and trying to create another (different) instance")
    {
        unsigned int length = 345;
        unsigned int q = 23;
        double p_mut = 0.6;
        double p_error = p_mut / 10;
        double p_effect = 0.2;
        double p_epistasis = 0.1;
        constants::Constants &cons = constants::Constants::create_instance(length, q, p_mut, p_error, p_effect,
                                                                           p_epistasis, std::filesystem::path());
        //        REQUIRE(cons.L == length);
        //        REQUIRE(cons.Q == q);
        //        REQUIRE(cons.P_MUT == p_mut);
        //        REQUIRE(cons.P_EFFECT == p_effect);
        //        REQUIRE(cons.P_EPISTASIS == p_epistasis);
        //        REQUIRE(cons.OUTPUT_DIR.empty());
        //
        //        length = 29348;
        //        q = 14;
        //        p_mut = 0.3;
        //        p_error = p_mut/10;
        //        p_effect = 0.7;
        //        p_epistasis = 0.1;
        //        constants::Constants &cons2 = constants::Constants::create_instance(length, q, p_mut, p_error,
        //        p_effect, p_epistasis, std::filesystem::path());
        REQUIRE(cons.L != length);
        REQUIRE(cons.Q != q);
        REQUIRE(cons.P_MUT != p_mut);
        REQUIRE(cons.P_EFFECT != p_effect);
        REQUIRE(cons.P_EPISTASIS != p_epistasis);

        REQUIRE(cons.OUTPUT_DIR.empty());
    }
}