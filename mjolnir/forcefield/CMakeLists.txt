set(mjolnir_forcefield_cpp_files)

add_subdirectory(3SPN2)
add_subdirectory(PDNS)
set(mjolnir_cpp_files ${mjolnir_cpp_files} ${mjolnir_forcefield_cpp_files} PARENT_SCOPE)
