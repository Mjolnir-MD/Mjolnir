# ----------------------------------------------------------------------------
# add executables

# message(STATUS ${mjolnir_source_files})

add_executable(mjolnir mjolnir.cpp)
set_target_properties(mjolnir PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
    COMPILE_FLAGS "${MJOLNIR_OPTIMIZATON_FLAGS}")

if(SEPARATE_BUILD)
    add_library(mjolnir_core STATIC ${mjolnir_source_files})
    set_target_properties(mjolnir_core PROPERTIES
        COMPILE_FLAGS "${MJOLNIR_OPTIMIZATON_FLAGS}")
    target_link_libraries(mjolnir mjolnir_core)
endif()

if(OpenMP_FOUND AND USE_OPENMP)
    message(STATUS "adding OpenMP flags ${OpenMP_CXX_FLAGS} to mjolnir ...")
    target_link_libraries(mjolnir ${OpenMP_CXX_LIBRARIES})

    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
        # After CMake 3.13, we can use target_link_options
        set_target_properties(mjolnir PROPERTIES LINK_FLAGS "-parallel")
    endif()
else()
    message(STATUS "Ignoring OpenMP ...")
endif()
