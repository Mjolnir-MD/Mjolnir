set(mjolnir_forcefield_local_cpp_files
    # -------------------------------------------------------------------------
    # potential
    "${CMAKE_CURRENT_SOURCE_DIR}/ClementiDihedralPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/CosinePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaussianPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GoContactPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GoContactAttractivePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GoContactRepulsivePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/HarmonicPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/PeriodicGaussianPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/WormLikeChainPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/WormLikeChainOffsetPotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/MBasinAttractivePotential.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/MBasinRepulsivePotential.cpp"
    # -------------------------------------------------------------------------
    # interaction
    "${CMAKE_CURRENT_SOURCE_DIR}/BondAngleInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BondLengthInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BondLengthGoContactInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DihedralAngleInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ContactInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GoContactInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DirectionalContactInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DummyInteraction.cpp"
    )
set(mjolnir_forcefield_cpp_files
    ${mjolnir_forcefield_cpp_files}
    ${mjolnir_forcefield_local_cpp_files} PARENT_SCOPE)
