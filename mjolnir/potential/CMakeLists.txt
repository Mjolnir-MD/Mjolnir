set(mjolnir_potential_cpp_files)

add_subdirectory(external)
add_subdirectory(global)
add_subdirectory(local)
set(mjolnir_cpp_files ${mjolnir_potential_cpp_files} PARENT_SCOPE)
