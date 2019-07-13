# ----------------------------------------------------------------------------
# add executables

# message(STATUS ${mjolnir_source_files})

add_executable(mjolnir mjolnir.cpp)
set_target_properties(mjolnir PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

if(SEPARATE_BUILD)
    add_library(mjolnir_core STATIC ${mjolnir_source_files})
    target_link_libraries(mjolnir mjolnir_core)
endif()

if(OpenMP_FOUND AND USE_OPENMP)
    message(STATUS "adding OpenMP flags ${OpenMP_CXX_FLAGS} to mjolnir ...")
    target_link_libraries(mjolnir ${OpenMP_CXX_LIBRARIES})
else()
    message(STATUS "Ignoring OpenMP ...")
endif()
