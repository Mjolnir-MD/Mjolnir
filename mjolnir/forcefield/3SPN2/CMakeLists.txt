set(mjolnir_forcefield_3SPN2_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/ThreeSPN2BasePairPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ThreeSPN2BasePairInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ThreeSPN2BaseStackingPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ThreeSPN2BaseStackingInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ThreeSPN2BondPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ThreeSPN2ExcludedVolumePotential.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_3SPN2_cpp_files} PARENT_SCOPE)
