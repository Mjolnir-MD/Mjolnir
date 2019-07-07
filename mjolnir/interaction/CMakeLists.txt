set(mjolnir_interaction_cpp_files)

# add_subdirectory(external)
# add_subdirectory(global)
add_subdirectory(local)
set(mjolnir_cpp_files ${mjolnir_cpp_files} ${mjolnir_interaction_cpp_files} PARENT_SCOPE)
