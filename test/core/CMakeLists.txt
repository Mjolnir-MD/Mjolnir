set(TEST_NAMES
    test_matrix
    test_vector
    test_math
    test_range
    test_throw_exception
    test_topology
    test_exclusion_list
    test_system_motion_remover
    test_file_inclusion
    test_fixed_vector
    test_parse_range

    test_harmonic_potential
    test_gaussian_potential
    test_periodic_gaussian_potential
    test_cosine_potential
    test_worm_like_chain_potential
    test_worm_like_chain_offset_potential
    test_flp_angle_potential
    test_flp_dihedral_potential
    test_clementi_dihedral_potential
    test_go_contact_potential
    test_go_attractive_potential
    test_go_repulsive_potential
    test_mbasin_attractive_potential
    test_mbasin_repulsive_potential
    test_3spn2_bond_potential
    test_3spn2_base_stacking_potential
    test_3spn2_base_pair_potential
    test_3spn2_excluded_volume_potential
    test_lennard_jones_potential
    test_uniform_lennard_jones_potential
    test_lennard_jones_attractive_potential
    test_wca_potential
    test_isolf_potential
    test_excluded_volume_potential
    test_inverse_power_potential
    test_hard_core_excluded_volume_potential
    test_debye_huckel_potential
    test_pdns_potential
    test_pwmcos_potential
    test_implicit_membrane_potential
    test_excluded_volume_wall_potential
    test_lennard_jones_wall_potential
    test_axis_aligned_plane_shape

    test_bond_length_interaction
    test_bond_length_gocontact_interaction
    test_contact_interaction
    test_contact_gocontact_interaction
    test_directional_contact_interaction
    test_bond_angle_interaction
    test_dihedral_angle_interaction
    test_3spn2_base_stacking_interaction
    test_3spn2_base_pair_interaction
    test_3spn2_cross_stacking_interaction
    test_global_pair_interaction
    test_global_pair_excluded_volume_interaction
    test_global_pair_lennard_jones_interaction
    test_pdns_interaction
    test_pwmcos_interaction
    test_external_distance_interaction
    test_position_restraint_interaction
    test_afm_fitting_interaction
    test_rectangular_box_interaction
    test_pulling_force_interaction

    test_multiple_basin_forcefield
    test_hybrid_forcefield

    test_neighbor_list
    test_unlimited_verlet_list
    test_periodic_verlet_list
    test_unlimited_cell_list
    test_periodic_cell_list
    test_zorder_rtree

    test_read_harmonic_potential
    test_read_go_contact_potential
    test_read_go_attractive_potential
    test_read_go_repulsive_potential
    test_read_mbasin_attractive_potential
    test_read_mbasin_repulsive_potential
    test_read_clementi_dihedral_potential
    test_read_gaussian_potential
    test_read_worm_like_chain_potential
    test_read_worm_like_chain_offset_potential
    test_read_periodic_gaussian_potential
    test_read_flexible_local_angle_potential
    test_read_flexible_local_dihedral_potential
    test_read_cosine_potential
    test_read_3spn2_bond_potential
    test_read_excluded_volume_potential
    test_read_inverse_power_potential
    test_read_hard_core_excluded_volume_potential
    test_read_3spn2_excluded_volume_potential
    test_read_debye_huckel_potential
    test_read_lennard_jones_potential
    test_read_uniform_lennard_jones_potential
    test_read_lennard_jones_attractive_potential
    test_read_tabulated_lennard_jones_attractive_potential
    test_read_wca_potential
    test_read_tabulated_wca_potential
    test_read_isolf_potential
    test_read_excluded_volume_wall_potential
    test_read_lennard_jones_wall_potential
    test_read_implicit_membrane_potential
    test_read_uniform_potential

    test_read_ignore_molecule
    test_read_ignore_group

    test_read_bond_length_interaction
    test_read_bond_angle_interaction
    test_read_dihedral_angle_interaction
    test_read_directional_contact_interaction
    test_read_dummy_interaction
    test_read_3spn2_base_stacking_interaction
    test_read_3spn2_base_pair_interaction
    test_read_global_pair_interaction
    test_read_pdns_interaction
    test_read_pwmcos_interaction
    test_read_position_restraint_interaction
    test_read_afm_fitting_interaction
    test_read_rectangular_box_interaction
    test_read_pulling_force_interaction
    test_read_com_pulling_force_interaction

    test_read_index_offset

    test_read_local_forcefield
    test_read_global_forcefield
    test_read_external_forcefield
    test_read_forcefield
    test_read_multiple_basin_forcefield
    test_read_hybrid_forcefield

    test_read_molecular_dynamics_simulator
    test_read_simulated_annealing_simulator
    test_read_steepest_descent_simulator
    test_read_switching_forcefield_simulator

    test_read_spatial_partition
    test_read_observer
    test_read_integrator
    test_read_system
    test_read_unit
    test_read_path

    test_save_load_msgpack
    )

if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message(STATUS "boost_unit_test_framework precompiled library found -> ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    # here we use less-aggressive optimization flag to turn on NaN checking.
    # -Ofast sometimes skips NaN check by replacing isnan always false.
    set_target_properties(${TEST_NAME} PROPERTIES
        COMPILE_FLAGS "${MJOLNIR_WARNING_FLAGS} -O2")

    if(SEPARATE_BUILD)
        target_link_libraries(${TEST_NAME} mjolnir_core)
    endif()
    if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
        target_link_libraries(${TEST_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    endif()

    if(OpenMP_CXX_FOUND AND USE_OPENMP)
        if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
            # After CMake 3.13, we can use target_link_options
            set_target_properties(${TEST_NAME} PROPERTIES LINK_FLAGS "-parallel")
        endif()
    endif()

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endforeach(TEST_NAME)
