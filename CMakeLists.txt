cmake_minimum_required (VERSION 2.6)

project(ABM)

set(CMAKE_C_FLAGS " -W -std=c++11 -O3")
set(CMAKE_CXX_FLAGS " -W -std=c++11 -O3")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
set(SOURCE_FILES src/main.cpp src/tinyxml2.cpp src/GrowthModel.h src/GrowthModel.cpp
                 src/Experiment.cpp src/Experiment.h src/ExperimentConstantGrowth.cpp
                 src/ExperimentConstantGrowth.h src/Cell.cpp src/Cell.h src/SimulationBuilder.cpp
                 src/SimulationBuilder.h src/ExperimentIteratedGrowth.cpp src/ExperimentIteratedGrowth.h)



find_package(Boost COMPONENTS system filesystem iostreams random REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})


add_executable(simulator ${SOURCE_FILES})
target_link_libraries( simulator ${Boost_LIBRARIES} )

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif(DOXYGEN_FOUND)
