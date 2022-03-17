set(mjolnir_forcefield_hybrid_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/HybridForceField.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DynamicHybridForceField.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_hybrid_cpp_files} PARENT_SCOPE)
