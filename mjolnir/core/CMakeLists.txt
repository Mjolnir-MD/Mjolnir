set(mjolnir_core_cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/System.cpp"
    )
set(mjolnir_cpp_files
    ${mjolnir_cpp_files}
    ${mjolnir_core_cpp_files} PARENT_SCOPE)
