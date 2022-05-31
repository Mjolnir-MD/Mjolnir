set(mjolnir_forcefield_stoichiometric_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalStoichiometricInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/StoichiometricCubicPanPotential.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_stoichiometric_cpp_files} PARENT_SCOPE)
