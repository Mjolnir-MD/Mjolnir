if(BUILD_UNIT_TEST)
    add_subdirectory(core)
endif()

if(BUILD_INTEGRATION_TEST)
    add_subdirectory(integration)
endif()

if(BUILD_UNIT_TEST AND USE_OPENMP AND OpenMP_CXX_FOUND)
    message(STATUS "building tests for OpenMP")
    add_subdirectory(omp)
endif()
