cmake_minimum_required(VERSION 3.0)
enable_testing()
project(mjolnir)

# the reference says it enables compiler's extension also (e.g. -std=gnu11,
# not -std=c++11). it's not so nice, so just add `std=c++11` instead.
# set(${CMAKE_CXX_STANDARD} 11)
# later we may need to prepare nicer stuff to set compilation flags depending
# on the environment.

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)

option(DEBUG "dump debug information" OFF)
if(DEBUG)
    add_definitions(-DMJOLNIR_DEBUG)
endif()

option(USE_APPROX "use approximation instruction for time-consuming math functions")
if(USE_APPROX)
    add_definitions(-DMJOLNIR_WITH_APPROX)
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extlib/toml/toml.hpp")
    execute_process(COMMAND git submodule update --init --recursive
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endif()

if(BOOST_ROOT)
    # if BOOST_ROOT is manually defined, use the library.
    find_package(Boost 1.67.0 REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    if(EXISTS "${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0/boost/version.hpp")
        message(STATUS "boost 1.69.0 exists.")
        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0)
    else()

        find_program(WGET wget)
        if(NOT WGET)
            message(FATAL_ERROR "wget is not installed. Can't download BOOST library.")
        endif()

        message(STATUS "downloading Boost to ./extlib/ ...")
        execute_process(COMMAND "${WGET}" https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.bz2
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            RESULT_VARIABLE DOWNLOAD_BOOST_RESULT OUTPUT_QUIET ERROR_QUIET)
        if(NOT ${DOWNLOAD_BOOST_RESULT} EQUAL "0")
            message(FATAL_ERROR "failed to download Boost. please confirm network connections.")
        endif()

        execute_process(COMMAND which sha256sum RESULT_VARIABLE SHA256SUM_EXISTS OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND which shasum    RESULT_VARIABLE SHASUM_EXISTS    OUTPUT_QUIET ERROR_QUIET)
        if(${SHA256SUM_EXISTS} EQUAL "0")
            execute_process(COMMAND sha256sum --check boost_1_69_0_tar_bz2_sha256sum.dat
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
                RESULT_VARIABLE   DOWNLOADED_SHA256SUM_RESULT)
        elseif(${SHASUM_EXISTS} EQUAL "0")
            execute_process(COMMAND shasum --algorithm 256 --check boost_1_69_0_tar_bz2_sha256sum.dat
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
                RESULT_VARIABLE   DOWNLOADED_SHA256SUM_RESULT)
        endif()

        if(NOT ${DOWNLOADED_SHA256SUM_RESULT} EQUAL "0")
            message(FATAL_ERROR "boost_1_69_0.tar.bz2 has invalid sha256sum!")
        endif()

        execute_process(COMMAND tar xf boost_1_69_0.tar.bz2
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            OUTPUT_QUIET ERROR_QUIET)

        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0)
        message(STATUS "done.")
    endif()
endif()

add_subdirectory(src)
add_subdirectory(test)
