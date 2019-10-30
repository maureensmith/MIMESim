

#TODO kann man das schöner mache mit den Pfaden?
set(DCABECHMARK_HEADER_LIST
        ../include/DCABenchmark/Constants.hpp
        ../include/DCABenchmark/Mutation.hpp
        ../include/DCABenchmark/Species.hpp
        ../include/DCABenchmark/Utils.hpp
        ../include/DCABenchmark/FunctionalSequence.hpp
        ../include/DCABenchmark/BindingCompetition.hpp
        )

set(DCABENCHMARK_SOURCE_LIST
        Constants.cpp
        Mutation.cpp
        Species.cpp
        Utils.cpp
        FunctionalSequence.cpp
        BindingCompetition.cpp
        )

#Set(BUILD_SHARED_LIBS ON)
#MESSAGE(HALLOOO ${BUILD_SHARED_LIBS})

#create a library for test purposes (a static one for now)
add_library(DCABenchmark ${DCABENCHMARK_SOURCE_LIST} ${DCABECHMARK_HEADER_LIST})
# tell the library that the headers are given in this directory
#PUBLICH: use it for specified target and for targets which links with this project
target_include_directories(DCABenchmark PUBLIC ../include/DCABenchmark)

# shared libraries need PIC
#set_property(TARGET DCABenchmark PROPERTY POSITION_INDEPENDENT_CODE 1)
target_link_libraries(DCABenchmark Eigen3::Eigen)
target_link_libraries(DCABenchmark stdc++fs)
target_link_libraries(DCABenchmark sam2counts)


add_executable(DCABenchmark_prog
        main.cpp
        #${DCABENCHMARK_SOURCE_LIST}
        #${DCABECHMARK_HEADER_LIST}
        )
# We need this directory, and users of our library will need it too
#target_include_directories(DCABenchmark_prog PUBLIC ../include)

target_link_libraries(DCABenchmark_prog DCABenchmark)


