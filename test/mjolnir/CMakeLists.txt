set(TEST_NAMES
    test_matrix
    test_vector
    test_range
    test_access_iterator
    test_throw_exception
    test_topology
    test_fast_inv_sqrt
    test_harmonic_potential
    test_gaussian_potential
    test_flp_angle_potential
    test_flp_dihedral_potential
    test_clementi_dihedral_potential
    test_go1012_potential
    test_lennard_jones_potential
    test_excluded_volume_potential
    test_debye_huckel_potential
    test_bond_length_interaction
    test_bond_angle_interaction
    test_dihedral_angle_interaction
    test_implicit_membrane_potential
    test_axis_aligned_plane_shape
    )

add_definitions("-O3 -std=c++11")

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
