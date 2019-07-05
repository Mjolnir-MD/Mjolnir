# ----------------------------------------------------------------------------
# add executables

add_executable(mjolnir mjolnir.cpp ${mjolnir_source_files})
set_target_properties(mjolnir PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

if(OpenMP_FOUND AND USE_OPENMP)
    message(STATUS "adding OpenMP flags ${OpenMP_CXX_FLAGS} to mjolnir ...")
    target_compile_options(mjolnir PUBLIC -DMJOLNIR_WITH_OPENMP ${OpenMP_CXX_FLAGS})
    target_link_libraries(mjolnir ${OpenMP_CXX_LIBRARIES})
else()
    message(STATUS "Ignoring OpenMP ...")
endif()
