set(mjolnir_omp_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/System.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ForceField.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/RandomNumberGenerator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BAOABLangevinIntegrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/UnderdampedLangevinIntegrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BondLengthInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ContactInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BondAngleInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DihedralAngleInteraction.cpp"
    )
set(mjolnir_cpp_files
    ${mjolnir_cpp_files}
    ${mjolnir_omp_cpp_files} PARENT_SCOPE)
