set(TEST_NAMES
    test_matrix
    test_vector
    test_range
    test_index_sequence
    test_zip_iterator
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
    test_zaxis_external_force_interaction
    )

add_definitions("-O3 -std=c++11")

find_library(BOOST_UNIT_TEST_FRAMEWORK_LIBRARY boost_unit_test_framework
             PATHS ${Boost_LIBRARY_DIRS} NO_DEFAULT_PATH)
if (BOOST_UNIT_TEST_FRAMEWORK_LIBRARY)
    message(STATUS "boost_unit_test_framework found @ ${BOOST_UNIT_TEST_FRAMEWORK_LIBRARY}")
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
else()
    message(STATUS "boost_unit_test_framework could not be found.")
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME} ${BOOST_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_NAME)
