add_subdirectory(core)
add_subdirectory(integration)

if(USE_OPENMP AND OpenMP_CXX_FOUND)
    message(STATUS "building tests for OpenMP")
    add_subdirectory(omp)
endif()
