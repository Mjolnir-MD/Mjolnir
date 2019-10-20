set(mjolnir_forcefield_PDNS_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/ProteinDNANonSpecificInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ProteinDNANonSpecificPotential.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_PDNS_cpp_files} PARENT_SCOPE)
