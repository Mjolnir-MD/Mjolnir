set(mjolnir_forcefield_PWMcos_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/PWMcosPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/PWMcosInteraction.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_PWMcos_cpp_files} PARENT_SCOPE)
