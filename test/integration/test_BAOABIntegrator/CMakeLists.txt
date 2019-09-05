# ----------------------------------------------------------------------------
# add executables

add_executable(test_BAOABIntegrator test_BAOABIntegrator.cpp)

set_target_properties(test_BAOABIntegrator PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

if(SEPARATE_BUILD)
    target_link_libraries(test_BAOABIntegrator mjolnir_core)
endif()

if(OpenMP_FOUND AND USE_OPENMP)
    message(STATUS "adding OpenMP flags ${OpenMP_CXX_FLAGS} to test_BAOABIntegrator ...")
    target_link_libraries(test_BAOABIntegrator ${OpenMP_CXX_LIBRARIES})

    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        # After CMake 3.13, we can use target_link_options
        set_target_properties(test_BAOABIntegrator PROPERTIES LINK_FLAGS "-parallel")
    endif()
endif()
