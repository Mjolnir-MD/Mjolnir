add_subdirectory(core)

if(USE_OPENMP AND OpenMP_CXX_FOUND)
    message(STATUS "building tests for OpenMP")
    add_subdirectory(omp)
endif()
