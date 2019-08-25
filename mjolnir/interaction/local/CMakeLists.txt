set(mjolnir_interaction_local_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/BondAngleInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BondLengthInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BondLengthGoContactInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DihedralAngleInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ContactInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GoContactInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DummyInteraction.cpp"
    )
set(mjolnir_interaction_cpp_files
    ${mjolnir_interaction_cpp_files}
    ${mjolnir_interaction_local_cpp_files} PARENT_SCOPE)
