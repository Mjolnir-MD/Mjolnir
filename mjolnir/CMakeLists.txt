set(mjolnir_cpp_files)

add_subdirectory(core)
add_subdirectory(input)
add_subdirectory(interaction)
add_subdirectory(potential)
if(USE_OPENMP AND OpenMP_CXX_FOUND)
    add_subdirectory(omp)
endif()

set(mjolnir_source_files ${mjolnir_source_files} ${mjolnir_cpp_files} PARENT_SCOPE)
