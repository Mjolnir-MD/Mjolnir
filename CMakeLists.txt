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

check_include_file_cxx(xmmintrin.h FOUND_SSE_INTRINSIC)
if(FOUND_SSE_INTRINSIC)
    check_cxx_source_compiles(
    "#include <xmmintrin.h>
     int main(){
         float r=1.;
         _mm_store_ss(&r, _mm_rsqrt_ss(_mm_load_ss(&r)));
         return 0;
     }" SSE_COMPILE_SUCCESS)

    if(SSE_COMPILE_SUCCESS)
        add_definitions(-DMJOLNIR_HAVE_SSE)
        add_definitions(-mmmx -msse)
    endif(SSE_COMPILE_SUCCESS)
endif(FOUND_SSE_INTRINSIC)

check_include_file_cxx(immintrin.h FOUND_AVX_INTRINSIC)
if(FOUND_AVX_INTRINSIC)
    set(CMAKE_REQUIRED_DEFINITIONS -mavx)
    check_cxx_source_compiles(
    "#include <immintrin.h>
     int main(){
         double r=1.;
         double xs[4];
         _mm256_storeu_pd(xs, _mm256_add_pd(_mm256_set1_pd(r), _mm256_set1_pd(r)));
         return 0;
     }" AVX_COMPILE_SUCCESS)

    if(AVX_COMPILE_SUCCESS)
        add_definitions(-DMJOLNIR_HAVE_AVX)
        add_definitions(-mavx)
    endif(AVX_COMPILE_SUCCESS)

    set(CMAKE_REQUIRED_DEFINITIONS -mavx2)
    check_cxx_source_compiles(
    "#include <immintrin.h>
     int main(){
         __m256i x = _mm256_add_epi32(_mm256_set1_epi32(1), _mm256_set1_epi32(2));
         return 0;
     }" AVX2_COMPILE_SUCCESS)

    if(AVX2_COMPILE_SUCCESS)
        add_definitions(-DMJOLNIR_HAVE_AVX2)
        add_definitions(-mavx2)
    endif(AVX2_COMPILE_SUCCESS)
endif(FOUND_AVX_INTRINSIC)


add_subdirectory(src)
add_subdirectory(test)
