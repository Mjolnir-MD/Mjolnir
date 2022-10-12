set(mjolnir_forcefield_external_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/ExcludedVolumeWallPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ImplicitMembranePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/LennardJonesWallPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/HarmonicGroovePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/RectangularBoxInteraction.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_external_cpp_files} PARENT_SCOPE)
