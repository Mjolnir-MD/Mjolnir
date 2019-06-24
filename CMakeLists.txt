cmake_minimum_required(VERSION 3.9) # FindOpenMP
enable_testing()
project(mjolnir)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD          11) # use c++11
set(CMAKE_CXX_EXTENSIONS       OFF) # not gnu11, c++11

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)

# -----------------------------------------------------------------------------
# find compilation flags

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

# -----------------------------------------------------------------------------
# check input

option(DEBUG "dump debug information" OFF)
if(DEBUG)
    add_definitions(-DMJOLNIR_DEBUG)
endif()

option(USE_APPROX "use approximation instruction for time-consuming math functions" OFF)
if(USE_APPROX)
    add_definitions(-DMJOLNIR_WITH_APPROX)
endif()

option(USE_OPENMP "use OpenMP for CPU parallelization" ON)
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP ${OpenMP_CXX_VERSION} found")
    else()
        message(WARNING "OpenMP not found. OpenMP support is disabled.")
    endif()
endif()

# -----------------------------------------------------------------------------

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
        message(STATUS "Boost library not found. Donwloading...")

        file(DOWNLOAD https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.bz2
            ${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0.tar.bz2
            EXPECTED_HASH SHA256=8f32d4617390d1c2d16f26a27ab60d97807b35440d45891fa340fc2648b04406
            STATUS MJOLNIR_DOWNLOAD_BOOST_STATUS)

        # check status
        list(GET MJOLNIR_DOWNLOAD_BOOST_STATUS 0 MJOLNIR_DOWNLOAD_BOOST_RESULT)
        if(NOT ${MJOLNIR_DOWNLOAD_BOOST_RESULT} EQUAL 0)
            list(GET MJOLNIR_DOWNLOAD_BOOST_STATUS 1 MJOLNIR_DOWNLOAD_BOOST_ERROR)
            message(FATAL_ERROR "failed to download Boost library."
                ${MJOLNIR_DOWNLOAD_BOOST_ERROR})
        endif()

        message(STATUS "downloading completed. Unpacking...")

        execute_process(COMMAND tar xf boost_1_69_0.tar.bz2
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
            OUTPUT_QUIET ERROR_QUIET)

        include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_69_0)
        message(STATUS "done.")
    endif()
endif()

add_subdirectory(src)
add_subdirectory(test)
