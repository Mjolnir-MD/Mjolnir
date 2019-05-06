add_executable(mjolnir mjolnir.cpp)
set_target_properties(mjolnir
    PROPERTIES
    COMPILE_FLAGS "-Ofast"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)

if(OpenMP_FOUND)
    target_link_libraries(mjolnir ${OpenMP_CXX_LIBRARIES})
endif()

# XXX we may later need to add a script to generate optimal compiler flags
#     depending on the environment.
