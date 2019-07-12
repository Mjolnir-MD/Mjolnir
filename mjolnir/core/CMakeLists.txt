set(mjolnir_core_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/System.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/RandomNumberGenerator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BAOABLangevinIntegrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/UnderdampedLangevinIntegrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/VelocityVerletIntegrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/MolecularDynamicsSimulator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SimulatedAnnealingSimulator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SteepestDescentSimulator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/UnlimitedGridCellList.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGridCellList.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/VerletList.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/NaivePairCalculation.cpp"
    )
set(mjolnir_cpp_files
    ${mjolnir_cpp_files}
    ${mjolnir_core_cpp_files} PARENT_SCOPE)
