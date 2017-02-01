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


include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)

check_include_file_cxx(xmmintrin.h FOUND_SSE)
if(FOUND_SSE)
    check_cxx_source_compiles(
    "#include <xmmintrin.h>
     int main(){
         float r=1.;
         _mm_store_ss(&r, _mm_rsqrt_ss(_mm_load_ss(&r)));
         return 0;
     }" SSE_COMPILE_SUCCESS)

    if(SSE_COMPILE_SUCCESS)
        add_definitions(-DHAVE_SSE)
        add_definitions(-mmmx -msse)
    endif(SSE_COMPILE_SUCCESS)
endif(FOUND_SSE)

add_subdirectory(src)
add_subdirectory(test)
