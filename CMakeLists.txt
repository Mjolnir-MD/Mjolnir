cmake_minimum_required(VERSION 3.0)
enable_testing()
project(mjolnir)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD          11) # use c++11
set(CMAKE_CXX_EXTENSIONS       OFF) # not gnu11, c++11

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)

option(DEBUG "dump debug information" OFF)
if(DEBUG)
    add_definitions(-DMJOLNIR_DEBUG)
endif()
option(USE_APPROX "use approximation instruction for time-consuming math functions" OFF)
if(USE_APPROX)
    add_definitions(-DMJOLNIR_WITH_APPROX)
endif()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wall"      MJOLNIR_COMPILER_HAS_WALL)
if(MJOLNIR_COMPILER_HAS_WALL)
    add_definitions(-Wall)
endif()
check_cxx_compiler_flag("-Wextra"    MJOLNIR_COMPILER_HAS_WEXTRA)
if(MJOLNIR_COMPILER_HAS_WEXTRA)
    add_definitions(-Wextra)
endif()
check_cxx_compiler_flag("-Wpedantic" MJOLNIR_COMPILER_HAS_WPEDANTIC)
if(MJOLNIR_COMPILER_HAS_WPEDANTIC)
    add_definitions(-Wpedantic)
endif()
check_cxx_compiler_flag("-march=native" MJOLNIR_COMPILER_HAS_MARCH_NATIVE)
if(MJOLNIR_COMPILER_HAS_MARCH_NATIVE)
    add_definitions(-march=native)
endif()
check_cxx_compiler_flag("-mtune=native" MJOLNIR_COMPILER_HAS_MTUNE_NATIVE)
if(MJOLNIR_COMPILER_HAS_MTUNE_NATIVE)
    add_definitions(-mtune=native)
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extlib/toml/toml.hpp")
    execute_process(COMMAND git submodule update --init --recursive
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endif()

option(FIND_BOOST "find boost library" OFF)
if(BOOST_ROOT OR FIND_BOOST)
    # if BOOST_ROOT is manually defined, use the library.
    find_package(Boost 1.67.0 REQUIRED COMPONENTS unit_test_framework)
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
