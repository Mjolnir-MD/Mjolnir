set(mjolnir_forcefield_MultipleBasin_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/MultipleBasin2BasinUnit.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/MultipleBasin3BasinUnit.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/MultipleBasinForceField.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/MBasinAttractivePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/MBasinRepulsivePotential.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_MultipleBasin_cpp_files} PARENT_SCOPE)
