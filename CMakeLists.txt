cmake_minimum_required(VERSION 3.9) # FindOpenMP
enable_testing()
project(mjolnir)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD          11) # use c++11
set(CMAKE_CXX_EXTENSIONS       OFF) # not gnu11, c++11

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)
add_definitions(-DTOML11_COLORIZE_ERROR_MESSAGE)

# =============================================================================
# check git hash and tag. then defime a macro, `MJOLNIR_VERSION_INFO`.

set(MJOLNIR_MOST_RECENT_TAG     "v0.0.0")
set(MJOLNIR_DIFFERENCE_FROM_TAG "0")
set(MJOLNIR_CURRENT_COMMIT_HASH "0")

if(APPLE)
    execute_process(
        COMMAND bash "-c" "git describe --tags --long HEAD | sed -E 's/(v[0-9]+\\.[0-9]+\\.[0-9]+)-([0-9]+)-g([a-z0-9]+)/\\1/'"
        OUTPUT_VARIABLE MJOLNIR_MOST_RECENT_TAG
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
    execute_process(
        COMMAND bash "-c" "git describe --tags --long HEAD | sed -E 's/(v[0-9]+\\.[0-9]+\\.[0-9]+)-([0-9]+)-g([a-z0-9]+)/\\2/'"
        OUTPUT_VARIABLE MJOLNIR_DIFFERENCE_FROM_TAG
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
    execute_process(
        COMMAND bash "-c" "git describe --tags --long HEAD | sed -E 's/(v[0-9]+\\.[0-9]+\\.[0-9]+)-([0-9]+)-g([a-z0-9]+)/\\3/'"
        OUTPUT_VARIABLE MJOLNIR_CURRENT_COMMIT_HASH
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
else()
    execute_process(
        COMMAND bash "-c" "git describe --tags --long HEAD | sed -r 's/(v[0-9]+\\.[0-9]+\\.[0-9]+)-([0-9]+)-g([a-z0-9]+)/\\1/'"
        OUTPUT_VARIABLE MJOLNIR_MOST_RECENT_TAG
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
    execute_process(
        COMMAND bash "-c" "git describe --tags --long HEAD | sed -r 's/(v[0-9]+\\.[0-9]+\\.[0-9]+)-([0-9]+)-g([a-z0-9]+)/\\2/'"
        OUTPUT_VARIABLE MJOLNIR_DIFFERENCE_FROM_TAG
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
    execute_process(
        COMMAND bash "-c" "git describe --tags --long HEAD | sed -r 's/(v[0-9]+\\.[0-9]+\\.[0-9]+)-([0-9]+)-g([a-z0-9]+)/\\3/'"
        OUTPUT_VARIABLE MJOLNIR_CURRENT_COMMIT_HASH
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_STRIP_TRAILING_WHITESPACE)
endif()

set(MJOLNIR_VERSION_INFO "${MJOLNIR_MOST_RECENT_TAG}")
if(NOT MJOLNIR_DIFFERENCE_FROM_TAG STREQUAL "0")
    set(MJOLNIR_VERSION_INFO "${MJOLNIR_MOST_RECENT_TAG}-dev (${MJOLNIR_CURRENT_COMMIT_HASH})")
endif()

add_definitions("-DMJOLNIR_VERSION=\"${MJOLNIR_VERSION_INFO}\"")
add_definitions("-DMJOLNIR_COMPILER_VERSION=\"${CMAKE_CXX_COMPILER}\"")

message(STATUS "MJOLNIR_VERSION = ${MJOLNIR_VERSION_INFO}")

# =============================================================================
# add warning flags

include(CheckCXXCompilerFlag)

set(MJOLNIR_WARNING_FLAGS "")

check_cxx_compiler_flag("-Wall"      MJOLNIR_COMPILER_HAS_WALL)
if(MJOLNIR_COMPILER_HAS_WALL)
    set(MJOLNIR_WARNING_FLAGS "${MJOLNIR_WARNING_FLAGS} -Wall")
endif()
check_cxx_compiler_flag("-Wextra"    MJOLNIR_COMPILER_HAS_WEXTRA)
if(MJOLNIR_COMPILER_HAS_WEXTRA)
    set(MJOLNIR_WARNING_FLAGS "${MJOLNIR_WARNING_FLAGS} -Wextra")
endif()
check_cxx_compiler_flag("-Wpedantic" MJOLNIR_COMPILER_HAS_WPEDANTIC)
if(MJOLNIR_COMPILER_HAS_WPEDANTIC)
    set(MJOLNIR_WARNING_FLAGS "${MJOLNIR_WARNING_FLAGS} -Wpedantic")
endif()

message(STATUS "warning flags are ${MJOLNIR_WARNING_FLAGS}")

# =============================================================================
# add optimization flags

set(MJOLNIR_OPTIMIZATION_FLAGS "")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    check_cxx_compiler_flag("-fast" MJOLNIR_COMPILER_HAS_FAST)
    if(MJOLNIR_COMPILER_HAS_FAST)
        set(MJOLNIR_OPTIMIZATION_FLAGS "${MJOLNIR_OPTIMIZATION_FLAGS} -fast")
    endif()
    check_cxx_compiler_flag("-xHost" MJOLNIR_COMPILER_HAS_XHOST)
    if(MJOLNIR_COMPILER_HAS_XHOST)
        set(MJOLNIR_OPTIMIZATION_FLAGS "${MJOLNIR_OPTIMIZATION_FLAGS} -xHost")
    endif()
    check_cxx_compiler_flag("-static-intel" MJOLNIR_COMPILER_HAS_STATIC_INTEL)
    if(MJOLNIR_COMPILER_HAS_STATIC_INTEL)
        set(MJOLNIR_OPTIMIZATION_FLAGS "${MJOLNIR_OPTIMIZATION_FLAGS} -static-intel")
    endif()
else()
    check_cxx_compiler_flag("-march=native" MJOLNIR_COMPILER_HAS_MARCH_NATIVE)
    if(MJOLNIR_COMPILER_HAS_MARCH_NATIVE)
        set(MJOLNIR_OPTIMIZATION_FLAGS "${MJOLNIR_OPTIMIZATION_FLAGS} -march=native")
    endif()
    check_cxx_compiler_flag("-mtune=native" MJOLNIR_COMPILER_HAS_MTUNE_NATIVE)
    if(MJOLNIR_COMPILER_HAS_MTUNE_NATIVE)
        set(MJOLNIR_OPTIMIZATION_FLAGS "${MJOLNIR_OPTIMIZATION_FLAGS} -mtune=native")
    endif()
    check_cxx_compiler_flag("-Ofast" MJOLNIR_COMPILER_HAS_OFAST)
    check_cxx_compiler_flag("-O3"    MJOLNIR_COMPILER_HAS_O3)
    if(MJOLNIR_COMPILER_HAS_OFAST)
        set(MJOLNIR_OPTIMIZATION_FLAGS "${MJOLNIR_OPTIMIZATION_FLAGS} -Ofast")
    elseif(MJOLNIR_COMPILER_HAS_O3)
        set(MJOLNIR_OPTIMIZATION_FLAGS "${MJOLNIR_OPTIMIZATION_FLAGS} -O3")
    endif()
endif()

message(STATUS "optimization flags are ${MJOLNIR_OPTIMIZATION_FLAGS}")

# =============================================================================
# check input flags

# -----------------------------------------------------------------------------
# Check debug flag

option(DEBUG "dump debug information" OFF)
if(DEBUG)
    add_definitions(-DMJOLNIR_DEBUG)
endif()

# -----------------------------------------------------------------------------
# check float/double flags

option(SINGLE_PRECISION_SUPPORT "support float (32-bit floating point) type" ON)
if(NOT SINGLE_PRECISION_SUPPORT)
    add_definitions(-DMJOLNIR_WITHOUT_SINGLE_PRECISION)
endif()

option(DOUBLE_PRECISION_SUPPORT "support double (64-bit floating point) type" ON)
if(NOT DOUBLE_PRECISION_SUPPORT)
    add_definitions(-DMJOLNIR_WITHOUT_DOUBLE_PRECISION)
endif()

if((NOT SINGLE_PRECISION_SUPPORT) AND (NOT DOUBLE_PRECISION_SUPPORT))
    message(FATAL_ERROR "Both float and double are disabled, nothing to be built")
endif()

# -----------------------------------------------------------------------------
# check boundary flags

option(UNLIMITED_BOUNDARY_SUPPORT "support unlimited boundary (infinite space)" ON)
if(NOT UNLIMITED_BOUNDARY_SUPPORT)
    add_definitions(-DMJOLNIR_WITHOUT_UNLIMITED_BOUNDARY)
endif()

option(PERIODIC_BOUNDARY_SUPPORT "support periodic boundary (cuboidal)" ON)
if(NOT PERIODIC_BOUNDARY_SUPPORT)
    add_definitions(-DMJOLNIR_WITHOUT_PERIODIC_BOUNDARY)
endif()

if((NOT UNLIMITED_BOUNDARY_SUPPORT) AND (NOT PERIODIC_BOUNDARY_SUPPORT))
    message(FATAL_ERROR "Both unlimited and periodic boundaries are disabled, nothing to be built")
endif()

# -----------------------------------------------------------------------------
# check approximation flag

option(USE_APPROX "\
use approximation instruction for time-consuming math functions" OFF)
if(USE_APPROX)
    add_definitions(-DMJOLNIR_WITH_APPROX)
endif()

# -----------------------------------------------------------------------------
# check openmp flag

option(USE_OPENMP "use OpenMP for CPU parallelization" ON)
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP ${OpenMP_CXX_VERSION} found")
        add_definitions(-DMJOLNIR_WITH_OPENMP)
        set(MJOLNIR_OPTIMIZATION_FLAGS "${MJOLNIR_OPTIMIZATION_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
        message(WARNING "OpenMP not found. OpenMP support is disabled.")
    endif()
endif()

# -----------------------------------------------------------------------------
# check whether unit/integration tests are needed (used later)

option(BUILD_INTEGRATION_TEST "\
build integration test that does not run automatically \
because it takes relatively long time" ON)

option(BUILD_UNIT_TEST "build unit tests" ON)

# -----------------------------------------------------------------------------
# check separate build

option(SEPARATE_BUILD "separate build process. It reduces the compilation time \
after building it once. But it possibly gains performance loss. So this is not \
for users, but for developers who compile this application many time with small\
 changes." OFF)

if(APPLE)
    # With apple clang, `~destructor() = defualt` does not work and causes a
    # linker error. So on apple, SEPARATE_BUILD option cannot be used.
    set(SEPARATE_BUILD OFF)
endif()

if(SEPARATE_BUILD)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        message(WARNING
            "-DSEPARATE_BUILD=ON option may slow down the runtime efficiency.\n"
            "This flag is for developers. If you don't recompile many times, "
            "remove the flag, delete the CMakeCache.txt and run cmake again.")
    endif()
    add_definitions(-DMJOLNIR_SEPARATE_BUILD)
#     add_definitions(-flto)
endif()

# =============================================================================
# download dependencies

execute_process(COMMAND git submodule update --init --recursive
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

option(FIND_BOOST "find boost library" OFF)
if(BUILD_UNIT_TEST)
    if(BOOST_ROOT OR FIND_BOOST)
        # if BOOST_ROOT is manually defined, use the library.
        find_package(Boost 1.67.0 REQUIRED COMPONENTS unit_test_framework)
        include_directories(${Boost_INCLUDE_DIRS})
    else()
        if(EXISTS "${PROJECT_SOURCE_DIR}/extlib/boost_1_73_0/boost/version.hpp")
            message(STATUS "boost 1.73.0 exists.")
            include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_73_0)
        else()
            message(STATUS "Boost library not found. Donwloading...")
            message(STATUS "It can take minutes. hang tight!")

            file(DOWNLOAD https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2
                ${PROJECT_SOURCE_DIR}/extlib/boost_1_73_0.tar.bz2
                EXPECTED_HASH SHA256=4eb3b8d442b426dc35346235c8733b5ae35ba431690e38c6a8263dce9fcbb402
                STATUS MJOLNIR_DOWNLOAD_BOOST_STATUS)

            # check status
            list(GET MJOLNIR_DOWNLOAD_BOOST_STATUS 0 MJOLNIR_DOWNLOAD_BOOST_RESULT)
            if(NOT ${MJOLNIR_DOWNLOAD_BOOST_RESULT} EQUAL 0)
                list(GET MJOLNIR_DOWNLOAD_BOOST_STATUS 1 MJOLNIR_DOWNLOAD_BOOST_ERROR)
                message(FATAL_ERROR "failed to download Boost library."
                    ${MJOLNIR_DOWNLOAD_BOOST_ERROR})
            endif()

            message(STATUS "downloading completed. Unpacking...")

            execute_process(COMMAND tar xf boost_1_73_0.tar.bz2
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/extlib"
                OUTPUT_QUIET ERROR_QUIET)

            include_directories(${PROJECT_SOURCE_DIR}/extlib/boost_1_73_0)
            message(STATUS "done.")
        endif()
    endif()
endif(BUILD_UNIT_TEST)

# =============================================================================
# Build the code

set(mjolnir_source_files)
if(SEPARATE_BUILD)
    add_subdirectory(mjolnir)
endif()

add_subdirectory(src)

if(BUILD_UNIT_TEST OR BUILD_INTEGRATION_TEST)
    add_subdirectory(test)
endif()
