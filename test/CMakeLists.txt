add_subdirectory(core)

if(BUILD_INTEGRATION_TEST)
    add_subdirectory(integration)
endif()

if(USE_OPENMP AND OpenMP_CXX_FOUND)
    message(STATUS "building tests for OpenMP")
    add_subdirectory(omp)
endif()
