add_executable(mjolnir mjolnir.cpp)
set_target_properties(mjolnir
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)
