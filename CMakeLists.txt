cmake_minimum_required(VERSION 3.10)
project(DCABenchmark)
#TODO Projekt umbeennnen: Ist ja erstmal nur Datenerzeugen, bzw Sequenz Evolution für gewisse Funktion/Evolutionary constraints simulieren:
set(DCABenchmark_VERSION_MAJOR 0)
set(DCABenchmark_VERSION_MINOR 1)
set(DCABenchmark_VERSION_PATCH 0)
set(DCABenchmark_VERSION "${DCABenchmark_VERSION_MAJOR}.${DCABenchmark_VERSION_MINOR}.${DCABenchmark_VERSION_PATCH}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${MAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")


#include_directories(extern DCABenchmark TestDCABenchmark)
#TODO ggf checken ob beim MIME + DCA überhaupt benötigt wird, dann hier rein schieben
include_directories(../extern)
#including static library for sam2counts
link_directories(../libs)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

#TODO nötig? Oder setzte ich das am Anfang im Projekt?
#set_property(TARGET DCABenchmark PROPERTY VERSION ${DCABenchmark_VERSION})

# Some other library that we import that was also built using CMake
# and has an exported target.
#find_package(sam2counts REQUIRED)


# The compiled library code is here
add_subdirectory(src)




################## Testing #########################################################

#if someone else includes this package, and they use BUILD_TESTING, they probably do not want the tests here
#if(CMAKE_PROJECT_NAME STREQUAL DCABenchmark AND BUILD_TESTING)
add_subdirectory(test)

    ##### TODO weg damit
#    target_include_directories(gbBase_Test PUBLIC ${GB_BASE_TEST_DIR})
#    target_compile_definitions(gbBase_Test PRIVATE $<$<CONFIG:Debug>:GHULBUS_CONFIG_ASSERT_LEVEL_DEBUG>)
#    target_compile_definitions(gbBase_Test PUBLIC $<$<CXX_COMPILER_ID:MSVC>:_ENABLE_EXTENDED_ALIGNED_STORAGE>)
#    add_test(NAME TestBase COMMAND gbBase_Test)
#endif()