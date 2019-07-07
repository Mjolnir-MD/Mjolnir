set(mjolnir_interaction_global_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairExcludedVolumeInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairLennardJonesInteraction.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GlobalPairUniformLennardJonesInteraction.cpp"
    )
set(mjolnir_interaction_cpp_files
    ${mjolnir_interaction_cpp_files}
    ${mjolnir_interaction_global_cpp_files} PARENT_SCOPE)
