set(TEST_NAMES
    test_matrix
    test_vector
    test_math
    test_unit
    test_range
    test_throw_exception
    test_static_string
    test_topology

    test_harmonic_potential
    test_gaussian_potential
    test_angular_gaussian_potential
    test_flp_angle_potential
    test_flp_dihedral_potential
    test_clementi_dihedral_potential
    test_go1012_potential
    test_lennard_jones_potential
    test_uniform_lennard_jones_potential
    test_excluded_volume_potential
    test_debye_huckel_potential
    test_bond_length_interaction
    test_bond_angle_interaction
    test_dihedral_angle_interaction
    test_implicit_membrane_potential
    test_axis_aligned_plane_shape

    test_read_harmonic_potential
    test_read_go1012_contact_potential
    test_read_clementi_dihedral_potential
    test_read_gaussian_potential
    test_read_angular_gaussian_potential
    test_read_flexible_local_angle_potential
    test_read_flexible_local_dihedral_potential
    test_read_excluded_volume_potential
    test_read_debye_huckel_potential
    test_read_lennard_jones_potential
    test_read_uniform_lennard_jones_potential
    test_read_excluded_volume_wall_potential
    test_read_lennard_jones_wall_potential
    test_read_implicit_membrane_potential

    test_read_bond_length_interaction
    test_read_bond_angle_interaction
    test_read_dihedral_angle_interaction
    test_read_global_pair_interaction

    )

add_definitions("-O2 -march=native -mtune=native -std=c++11")

if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message(STATUS "boost_unit_test_framework precompiled library found -> ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
        target_link_libraries(${TEST_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    endif (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_NAME)
