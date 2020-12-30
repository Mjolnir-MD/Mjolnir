# ----------------------------------------------------------------------------
# add executables

add_executable(test_G-JFIntegrator test_G-JFIntegrator.cpp)

set_target_properties(test_G-JFIntegrator PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

if(SEPARATE_BUILD)
    target_link_libraries(test_G-JFIntegrator mjolnir_core)
endif()

if(OpenMP_CXX_FOUND AND USE_OPENMP)
    message(STATUS "adding OpenMP flags ${OpenMP_CXX_FLAGS} to test_G-JFIntegrator ...")
    target_link_libraries(test_G-JFIntegrator ${OpenMP_CXX_LIBRARIES})

    # here use the same options as src/mjolnir because it is an integration test.
    set_target_properties(test_G-JFIntegrator PROPERTIES
        COMPILE_FLAGS "${MJOLNIR_WARNING_FLAGS} ${MJOLNIR_OPTIMIZATON_FLAGS}")

    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        # After CMake 3.13, we can use target_link_options
        set_target_properties(test_G-JFIntegrator PROPERTIES LINK_FLAGS "-parallel")
    endif()
endif()
