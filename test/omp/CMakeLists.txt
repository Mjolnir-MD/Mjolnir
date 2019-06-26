set(TEST_NAMES
    test_omp_sort
    test_omp_random_number_generator

    test_omp_bond_length_interaction
    test_omp_contact_interaction
    test_omp_bond_angle_interaction
    test_omp_dihedral_angle_interaction

    test_omp_global_lennard_jones_interaction
    test_omp_global_uniform_lennard_jones_interaction
    test_omp_global_excluded_volume_interaction
    test_omp_global_debye_huckel_interaction

    test_omp_external_distance_interaction
    test_omp_position_restraint_interaction
    )

if(NOT (OpenMP_FOUND AND USE_OPENMP))
    message(FATAL_ERROR "Test codes for OpenMP implementation requires OpenMP library.")
endif()

if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message(STATUS "boost_unit_test_framework precompiled library found -> ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)

    target_link_libraries(${TEST_NAME} ${OpenMP_CXX_LIBRARIES})
    target_compile_options(${TEST_NAME} PUBLIC -DMJOLNIR_WITH_OPENMP ${OpenMP_CXX_FLAGS})

    if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
        target_link_libraries(${TEST_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    endif (Boost_UNIT_TEST_FRAMEWORK_FOUND)

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} --log_level=warning WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endforeach(TEST_NAME)
