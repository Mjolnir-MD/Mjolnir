set(mjolnir_forcefield_FLP_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/FlexibleLocalAnglePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/FlexibleLocalDihedralPotential.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_FLP_cpp_files} PARENT_SCOPE)
