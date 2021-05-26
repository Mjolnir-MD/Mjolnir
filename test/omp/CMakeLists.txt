set(TEST_NAMES
    test_omp_sort
    test_omp_random_number_generator
    test_omp_save_load_msgpack
    test_omp_system_motion_remover
    test_omp_multiple_basin_forcefield

    test_omp_bond_length_interaction
    test_omp_bond_length_gocontact_interaction
    test_omp_contact_interaction
    test_omp_gocontact_interaction
    test_omp_bond_angle_interaction
    test_omp_dihedral_angle_interaction
    test_omp_3spn2_base_base_interaction

    test_omp_global_lennard_jones_interaction
    test_omp_global_uniform_lennard_jones_interaction
    test_omp_global_excluded_volume_interaction
    test_omp_global_debye_huckel_interaction
    test_omp_global_pdns_interaction
    test_omp_global_pwmcos_interaction

    test_omp_external_distance_interaction
    test_omp_position_restraint_interaction
    )

if(NOT (OpenMP_CXX_FOUND AND USE_OPENMP))
    message(FATAL_ERROR "Test codes for OpenMP implementation requires OpenMP library.")
endif()

if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message(STATUS "boost_unit_test_framework precompiled library found -> ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)

    # here we use less-aggressive optimization flag to turn on NaN checking.
    set_target_properties(${TEST_NAME} PROPERTIES
        COMPILE_FLAGS "${MJOLNIR_WARNING_FLAGS} -O2 ${OpenMP_CXX_FLAGS}")

    if(SEPARATE_BUILD)
        target_link_libraries(${TEST_NAME} mjolnir_core)
    endif()

    target_link_libraries(${TEST_NAME} ${OpenMP_CXX_LIBRARIES})

    if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
        target_link_libraries(${TEST_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    endif (Boost_UNIT_TEST_FRAMEWORK_FOUND)

    if(OpenMP_FOUND AND USE_OPENMP)
        if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
            # After CMake 3.13, we can use target_link_options
            set_target_properties(${TEST_NAME} PROPERTIES LINK_FLAGS "-parallel")
        endif()
    endif()

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} --log_level=warning WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endforeach(TEST_NAME)
