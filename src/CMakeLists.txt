# ----------------------------------------------------------------------------
# check optimization flags

check_cxx_compiler_flag("-march=native" MJOLNIR_COMPILER_HAS_MARCH_NATIVE)
if(MJOLNIR_COMPILER_HAS_MARCH_NATIVE)
    add_definitions(-march=native)
endif()
check_cxx_compiler_flag("-mtune=native" MJOLNIR_COMPILER_HAS_MTUNE_NATIVE)
if(MJOLNIR_COMPILER_HAS_MTUNE_NATIVE)
    add_definitions(-mtune=native)
endif()
check_cxx_compiler_flag("-Ofast" MJOLNIR_COMPILER_HAS_OFAST)
check_cxx_compiler_flag("-O3"    MJOLNIR_COMPILER_HAS_O3)
if(MJOLNIR_COMPILER_HAS_OFAST)
    add_definitions(-Ofast)
elseif(MJOLNIR_COMPILER_HAS_O3)
    add_definitions(-O3)
endif()

# ----------------------------------------------------------------------------
# add executables

add_executable(mjolnir mjolnir.cpp)
set_target_properties(mjolnir PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

if(OpenMP_FOUND AND USE_OPENMP)
    message(STATUS "adding OpenMP flags ${OpenMP_CXX_FLAGS} to mjolnir ...")
    target_compile_options(mjolnir PUBLIC -DMJOLNIR_WITH_OPENMP ${OpenMP_CXX_FLAGS})
    target_link_libraries(mjolnir ${OpenMP_CXX_LIBRARIES})
else()
    message(STATUS "Ignoring OpenMP ...")
endif()
