set(mjolnir_forcefield_iSoLF_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/iSoLFAttractivePotential.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_iSoLF_cpp_files} PARENT_SCOPE)
