cmake_minimum_required(VERSION 3.0)
enable_testing()
project(mjolnir)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)

option(DEBUG "dump debug information" OFF)
if(DEBUG)
    add_definitions(-DMJOLNIR_DEBUG)
endif()
set(MJOLNIR_EXECUTABLE_LOCATION "${PROJECT_SOURCE_DIR}/bin")

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extlib/toml/toml.hpp")
    execute_process(COMMAND git submodule update --init --recursive
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endif()

# if manually defined
if(BOOST_ROOT)
    find_package(Boost 1.66.0 REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    if(EXISTS "${PROJECT_SOURCE_DIR}/extlib/boost_1_67_0/boost/version.hpp")
        message(STATUS "boost 1.67.0 exists.")
        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_67_0)
    else()
        message(STATUS "downloading Boost 1.67.0 to ./extlib/ ...")
        execute_process(COMMAND wget https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.bz2
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            RESULT_VARIABLE DOWNLOAD_BOOST_RESULT OUTPUT_QUIET ERROR_QUIET)
        if(NOT ${DOWNLOAD_BOOST_RESULT} EQUAL "0")
            message(FATAL_ERROR "failed to download Boost 1.67.0. please confirm network connections.")
        endif()

        execute_process(COMMAND which sha256sum RESULT_VARIABLE SHA256SUM_EXISTS OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND which shasum    RESULT_VARIABLE SHASUM_EXISTS    OUTPUT_QUIET ERROR_QUIET)
        if(${SHA256SUM_EXISTS} EQUAL "0")
            execute_process(COMMAND sha256sum --check boost_1_67_0_tar_bz2_sha256sum.dat
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
                RESULT_VARIABLE   DOWNLOADED_SHA256SUM_RESULT)
        elseif(${SHASUM_EXISTS} EQUAL "0")
            execute_process(COMMAND shasum --algorithm 256 --check boost_1_67_0_tar_bz2_sha256sum.dat
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
                RESULT_VARIABLE   DOWNLOADED_SHA256SUM_RESULT)
        endif()

        if(NOT ${DOWNLOADED_SHA256SUM_RESULT} EQUAL "0")
            message(FATAL_ERROR "boost_1_67_0.tar.bz2 has invalid sha256sum!")
        endif()

        execute_process(COMMAND tar xf boost_1_67_0.tar.bz2
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            OUTPUT_QUIET ERROR_QUIET)

        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_67_0)
        message(STATUS "done.")
    endif()
endif()

add_definitions("-std=c++11 -Ofast -march=native -mtune=native")

add_subdirectory(src)
add_subdirectory(test)
