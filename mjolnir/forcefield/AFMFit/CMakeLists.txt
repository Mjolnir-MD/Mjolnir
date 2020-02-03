set(mjolnir_forcefield_AFMFit_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/AFMFitInteraction.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_AFMFit_cpp_files} PARENT_SCOPE)
