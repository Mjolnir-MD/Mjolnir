set(TEST_NAMES
    test_matrix
    test_vector
    test_fast_inv_sqrt
    test_lennard_jones_potential
    test_excluded_volume_potential
    test_harmonic_potential
    test_clementi_dihedral_potential
    test_go1012_potential
    test_bond_length_interaction
    test_bond_angle_interaction
    test_dihedral_angle_interaction
    )

add_definitions("-O3 -std=c++11")

set(test_library_dependencies)
find_library(BOOST_UNITTEST_FRAMEWORK_LIBRARY boost_unit_test_framework)
if (BOOST_UNITTEST_FRAMEWORK_LIBRARY)
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
    set(test_library_dependencies boost_unit_test_framework)
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME} ${test_library_dependencies})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_NAME)
