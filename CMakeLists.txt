cmake_minimum_required(VERSION 2.8)
enable_testing()
project(mjolnir)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)

option(DEBUG "dump debug information" OFF)
if(DEBUG)
    add_definitions(-DMJOLNIR_DEBUG)
else(DEBUG)
    add_definitions(-DMJOLNIR_NO_DEBUG)
endif(DEBUG)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extlib/toml/toml.hpp")
    execute_process(COMMAND git submodule update --init --recursive
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endif()

add_definitions("-std=c++11 -Ofast -march=native -mtune=native")

find_package(Boost COMPONENTS unit_test_framework)
if(${Boost_FOUND})
    message(STATUS "Boost include dir = ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost library dir = ${Boost_LIBRARY_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
else()
    message(STATUS "Boost couldn't be found. some functionalities are disabled.")
    message(STATUS "If you have Boost somewhere cmake didn't look (e.g. /opt/boost-1.??.??), try")
    message(STATUS "`cmake path/to/CMakeLists.txt -DBOOST_ROOT=/path/to/boost`")
endif()

add_subdirectory(src)
add_subdirectory(test)
