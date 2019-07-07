set(mjolnir_input_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/read_units.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/read_system.cpp"
    )
set(mjolnir_cpp_files
    ${mjolnir_cpp_files}
    ${mjolnir_input_cpp_files} PARENT_SCOPE)
