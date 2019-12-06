set(mjolnir_potential_local_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/ClementiDihedralPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/CosinePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/FlexibleLocalAnglePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/FlexibleLocalDihedralPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaussianPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GoContactPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GoContactAttractivePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GoContactRepulsivePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/HarmonicPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGaussianPotential.cpp"
    )
set(mjolnir_potential_cpp_files
    ${mjolnir_potential_cpp_files}
    ${mjolnir_potential_local_cpp_files} PARENT_SCOPE)
