set(mjolnir_input_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/read_input_file.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_integrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_local_forcefield.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_local_interaction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_local_potential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_observer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_simulator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_spatial_partition.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_system.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_units.cpp"
    )
set(mjolnir_cpp_files
    ${mjolnir_cpp_files}
    ${mjolnir_input_cpp_files} PARENT_SCOPE)
