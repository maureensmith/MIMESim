
find_package(Catch2 REQUIRED)
#TODO was genau macht das?
enable_testing()
#enables testing and set a BUILD_TESTING option so users can turn testing on and off
include(CTest)
include(Catch)

set(DCABENCHMARK_TEST_SOURCES
        catch.hpp
        test_bindingCompetition.cpp
        test_functionalSequence.cpp
        test_main.cpp  test_species.cpp
        test_utils.cpp
        )
# target for the test executable
add_executable(TestDCABenchmark ${DCABENCHMARK_TEST_SOURCES})
# Headers of the library to be tested
target_include_directories(TestDCABenchmark PUBLIC ../include/DCABenchmark)

# linking to the main library, as well as the Catch2 testing library
target_link_libraries(TestDCABenchmark PRIVATE
        DCABenchmark
        Catch2::Catch2
        Eigen3::Eigen
        stdc++fs
        sam2counts
        )

# add
add_test(DCABenchmark TestDCABenchmark) # Command can be a target

