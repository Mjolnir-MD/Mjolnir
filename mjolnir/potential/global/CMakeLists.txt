set(mjolnir_potential_global_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/DebyeHuckelPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ExcludedVolumePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/InversePowerPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/HardCoreExcludedVolumePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/LennardJonesPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/UniformLennardJonesPotential.cpp"
    )
set(mjolnir_potential_cpp_files
    ${mjolnir_potential_cpp_files}
    ${mjolnir_potential_global_cpp_files} PARENT_SCOPE)
