
#enable CTest
enable_testing()
#find_package(Catch2 REQUIRED)
#enables testing and set a BUILD_TESTING option so users can turn testing on and off
include(CTest)
#include(Catch)

# the sources are necessary for all tests
set(CATCH_SOURCES catch.hpp test_main.cpp)

set(DCABENCHMARK_TEST_SOURCES
        #test_constants.cpp
        test_bindingCompetition.cpp
        test_functionalSequence.cpp
        test_species.cpp
        test_utils.cpp
        test_mutation.cpp)

# target for the test of the constants file
add_executable(TestConstants ${CATCH_SOURCES} test_constants.cpp)

# target for the test executable
#add_executable(TestDCABenchmark ${CATCH_SOURCES} ${DCABENCHMARK_TEST_SOURCES})
add_executable(TestSpecies ${CATCH_SOURCES} test_species.cpp)



# Headers of the library to be tested
target_include_directories(TestSpecies PUBLIC ../include/DCABenchmark)
target_include_directories(TestConstants PUBLIC ../include/DCABenchmark)

# linking to the main library, as well as the Catch2 testing library
#target_link_libraries(TestSpecies PRIVATE DCABenchmark Catch2::Catch2)
#target_link_libraries(TestConstants PRIVATE DCABenchmark Catch2::Catch2)
target_link_libraries(TestSpecies PRIVATE DCABenchmark)
target_link_libraries(TestConstants PRIVATE DCABenchmark)



#TODO not necessary.... doch machen, pro laufzeit/singleton ein add_test

# The struct containing the constant parameters is given as singleton pattern.
# However, this is not easy to test, since, once an object has buin created, it is not possible to create a new one (also not in the different test cases and section or by restarting in the main)
# Hence: Different test for different singletons, i.e. parameter settings

# this requires a functional immaculate constant class.... test this first.

#add_test(TestSingletonBehaviour TestConstants)


# add test for the rest with q = 2
#add_test(TestWithQ2 TestDCABenchmark huhu) # Command can be a target

