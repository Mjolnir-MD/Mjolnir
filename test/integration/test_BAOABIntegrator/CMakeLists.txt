add_executable(test_BAOABIntegrator test_BAOABIntegrator.cpp)

set_target_properties(test_BAOABIntegrator PROPERTIES
    COMPILE_FLAGS "-Ofast"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

if(USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(test_BAOABIntegrator ${OpenMP_CXX_LIBRARIES})
endif()
