set(mjolnir_potential_external_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/ExcludedVolumeWallPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ImplicitMembranePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/LennardJonesWallPotential.cpp"
    )
set(mjolnir_potential_cpp_files
    ${mjolnir_potential_cpp_files}
    ${mjolnir_potential_external_cpp_files} PARENT_SCOPE)
