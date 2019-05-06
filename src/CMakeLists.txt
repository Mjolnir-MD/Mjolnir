add_executable(mjolnir mjolnir.cpp)
set_target_properties(mjolnir
    PROPERTIES
    COMPILE_FLAGS "-Ofast -march=native -mtune=native"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)

# XXX we may later need to add a script to generate optimal compiler flags
#     depending on the environment.
