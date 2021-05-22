set(mjolnir_forcefield_global_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairExcludedVolumeInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairLennardJonesInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DebyeHuckelPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ExcludedVolumePotential.cpp"
#     "${CMAKE_CURRENT_SOURCE_DIR}/InversePowerPotential.cpp"
#     "${CMAKE_CURRENT_SOURCE_DIR}/HardCoreExcludedVolumePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/LennardJonesPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/LennardJonesAttractivePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/WCAPotential.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_global_cpp_files} PARENT_SCOPE)
