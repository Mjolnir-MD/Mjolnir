#define BOOST_TEST_MODULE "test_omp_3spn2_base_base_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/math/math.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/ThreeSPN2BaseBaseInteraction.hpp>

BOOST_AUTO_TEST_CASE(omp_ThreeSPN2BasePairIntearction_basepair)
{
    namespace test = mjolnir::test;

    constexpr double tol = 1e-4;
    mjolnir::LoggerManager::set_default_logger(
            "test_omp_3spn2_base_base_interaction.log");

    using seq_traits_type  = mjolnir::SimulatorTraits      <double, mjolnir::UnlimitedBoundary>;
    using omp_traits_type  = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = typename omp_traits_type::real_type;
    using coordinate_type  = typename omp_traits_type::coordinate_type;
    using boundary_type    = typename omp_traits_type::boundary_type;
    using topology_type    = mjolnir::Topology;

    using potential_type       = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;

    using omp_system_type      = mjolnir::System<omp_traits_type>;
    using omp_interaction_type = mjolnir::ThreeSPN2BaseBaseInteraction<omp_traits_type>;

    using seq_system_type      = mjolnir::System<seq_traits_type>;
    using seq_interaction_type = mjolnir::ThreeSPN2BaseBaseInteraction<seq_traits_type>;

    using omp_parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<omp_traits_type>;
    using seq_parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<seq_traits_type>;
    using omp_parameter_type      = typename omp_parameter_list_type::parameter_type;
    using seq_parameter_type      = typename seq_parameter_list_type::parameter_type;

    using omp_partition_type      = mjolnir::VerletList<omp_traits_type, potential_type>;
    using seq_partition_type      = mjolnir::VerletList<seq_traits_type, potential_type>;

    // those actually does not depends on traits.
    using base_kind            = typename omp_parameter_list_type::base_kind;
    using ignore_group_type    = typename omp_parameter_list_type::ignore_group_type;
    using ignore_molecule_type = typename omp_parameter_list_type::ignore_molecule_type;

    constexpr real_type pi = mjolnir::math::constants<real_type>::pi();

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal() / 1000.0) *
                          unit_type::avogadro_constant());
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal()) /
                            unit_type::avogadro_constant());
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom());

        phys_type::set_m_to_length(unit_type::m_to_angstrom());
        phys_type::set_length_to_m(unit_type::angstrom_to_m());

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom(), 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m(), 3));

        phys_type::set_length_unit(length);
    }

    // Here, this test case checks base pairing interaction only.
    //
    //  theta1    theta2
    //       |    |
    // Si o  v    v o Sj
    //     \-.   .-/
    //      o === o
    //     Bi  r  Bj
    //
    // Required test cases are the combination of the following conditions
    //
    // BasePairing:
    // rij:
    //  1. (r < r0)
    //  2. (r0 < r)
    // theta1:
    //  1. theta < pi/2K
    //  2. pi/2K < theta < pi/K
    //  3. pi/K  < theta
    // theta2:
    //  1. theta < pi/2K
    //  2. pi/2K < theta < pi/K
    //  3. pi/K  < theta

    std::mt19937 rng(123456789);

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        const std::vector<base_kind> bases{base_kind::A, base_kind::T};
        //  theta1    theta2
        //       |    |
        //  0 o  v    v o 2
        //     \-.   .-/
        //      o === o
        //     1       3
        //
        constexpr auto invalid = potential_type::invalid();

        omp_parameter_list_type omp_parameter_list(
            mjolnir::ThreeSPN2BaseBaseGlobalPotentialParameter<real_type>{}, {
                {1, omp_parameter_type{bases.at(0), 0, 0, invalid, invalid}},
                {3, omp_parameter_type{bases.at(1), 1, 2, invalid, invalid}}
            }, {}, ignore_molecule_type("Nothing"), ignore_group_type({})
        );
        seq_parameter_list_type seq_parameter_list(
            mjolnir::ThreeSPN2BaseBaseGlobalPotentialParameter<real_type>{}, {
                {1, seq_parameter_type{bases.at(0), 0, 0, invalid, invalid}},
                {3, seq_parameter_type{bases.at(1), 1, 2, invalid, invalid}}
            }, {}, ignore_molecule_type("Nothing"), ignore_group_type({})
        );

        omp_interaction_type omp_interaction(potential_type{},
            omp_parameter_list_type(omp_parameter_list),
            mjolnir::SpatialPartition<omp_traits_type, potential_type>(
                mjolnir::make_unique<omp_partition_type>()));

        seq_interaction_type seq_interaction(potential_type{},
            seq_parameter_list_type(seq_parameter_list),
            mjolnir::SpatialPartition<seq_traits_type, potential_type>(
                mjolnir::make_unique<seq_partition_type>()));

        omp_system_type omp_sys(4, boundary_type{});
        seq_system_type seq_sys(4, boundary_type{});

        test::clear_everything(omp_sys);
        test::clear_everything(seq_sys);

        topology_type topol(4);
        topol.add_connection(0, 1, "bond");
        topol.add_connection(2, 3, "bond");
        topol.construct_molecules();

        omp_interaction.initialize(omp_sys, topol);
        seq_interaction.initialize(seq_sys, topol);

        // paramter_lists in the interactions will be initialized via
        // interaction.initialize(), but parameter_list here will not.
        potential_type pot;
        omp_parameter_list.initialize(omp_sys, topol, pot);

        const auto bp_kind      = omp_parameter_list.bp_kind(bases.at(0), bases.at(1));
        const auto rbp_0        = omp_parameter_list.r0(bp_kind);
        const auto theta1_0     = omp_parameter_list.theta1_0(bp_kind);
        const auto theta2_0     = omp_parameter_list.theta2_0(bp_kind);
        const auto pi_over_K_BP = omp_parameter_list.pi_over_K_BP();

        for(const auto rbp  : {rbp_0 - 0.2, rbp_0 + 0.5})
        {
        for(const auto theta1 : {theta1_0 + pi_over_K_BP * 0.2,
                                 theta1_0 + pi_over_K_BP * 0.7,
                                 theta1_0 + pi_over_K_BP * 1.2,
                                 theta1_0 - pi_over_K_BP * 1.2,
                                 theta1_0 - pi_over_K_BP * 0.7,
                                 theta1_0 - pi_over_K_BP * 0.2})
        {
        for(const auto theta2 : {theta2_0 + pi_over_K_BP * 0.2,
                                 theta2_0 + pi_over_K_BP * 0.7,
                                 theta2_0 + pi_over_K_BP * 1.2,
                                 theta2_0 - pi_over_K_BP * 1.2,
                                 theta2_0 - pi_over_K_BP * 0.7,
                                 theta2_0 - pi_over_K_BP * 0.2})
        {
            BOOST_TEST_MESSAGE("===========================================");
            BOOST_TEST_MESSAGE("rbp    = " << rbp);
            BOOST_TEST_MESSAGE("theta1 = " << theta1);
            BOOST_TEST_MESSAGE("theta2 = " << theta2);

            //  theta1    theta2
            //       |    |
            //  0 o  v    v o 2
            //     \-.   .-/
            //      o === o
            //     1       3
            //
            omp_sys.position(1) = coordinate_type(0.0, 0.0, 0.0);
            omp_sys.position(3) = coordinate_type(rbp, 0.0, 0.0);

            omp_sys.position(0) = coordinate_type(std::cos(theta1),          std::sin(theta1),    0.0);
            omp_sys.position(2) = coordinate_type(std::cos(pi-theta2) + rbp, std::sin(pi-theta2), 0.0);

            test::apply_random_rotation(omp_sys, rng);

            // check configuration is okay
            {
                const auto v13 = omp_sys.position(3) - omp_sys.position(1);
                BOOST_TEST_REQUIRE(mjolnir::math::length(v13) == rbp,
                                   boost::test_tools::tolerance(1e-3));

                const auto v10 = omp_sys.position(0) - omp_sys.position(1);
                const auto cos_013 = mjolnir::math::dot_product(v13, v10) /
                    (mjolnir::math::length(v13) * mjolnir::math::length(v10));
                BOOST_TEST_REQUIRE(std::acos(cos_013) == theta1,
                                   boost::test_tools::tolerance(1e-3));

                const auto v23 = omp_sys.position(3) - omp_sys.position(2);
                const auto cos_132 = mjolnir::math::dot_product(v13, v23) /
                    (mjolnir::math::length(v13) * mjolnir::math::length(v23));
                BOOST_TEST_REQUIRE(std::acos(cos_132) == theta2,
                                   boost::test_tools::tolerance(1e-3));
            }

            test::apply_random_perturbation(omp_sys, rng, 0.01);

            for(std::size_t i=0; i<seq_sys.size(); ++i)
            {
                seq_sys.position(i) = omp_sys.position(i);
            }
            test::clear_force(omp_sys);
            test::clear_force(seq_sys);

            test::check_force_consistency              (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
            test::check_force_and_energy_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
            test::check_force_and_virial_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
            test::check_force_energy_virial_consistency(omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
            test::check_energy_consistency             (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);

        } // theta2
        } // theta1
        } // rbp
    }
}

BOOST_AUTO_TEST_CASE(ThreeSPN2CrossStackingIntearction_numerical_diff)
{
    namespace test = mjolnir::test;

    constexpr double tol = 1e-4;
    mjolnir::LoggerManager::set_default_logger(
            "test_omp_3spn2_base_base_interaction.log");

    using seq_traits_type  = mjolnir::SimulatorTraits      <double, mjolnir::UnlimitedBoundary>;
    using omp_traits_type  = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = typename omp_traits_type::real_type;
    using matrix33_type    = typename omp_traits_type::matrix33_type;
    using coordinate_type  = typename omp_traits_type::coordinate_type;
    using boundary_type    = typename omp_traits_type::boundary_type;
    using topology_type    = mjolnir::Topology;

    using potential_type       = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;

    using omp_system_type      = mjolnir::System<omp_traits_type>;
    using omp_interaction_type = mjolnir::ThreeSPN2BaseBaseInteraction<omp_traits_type>;

    using seq_system_type      = mjolnir::System<seq_traits_type>;
    using seq_interaction_type = mjolnir::ThreeSPN2BaseBaseInteraction<seq_traits_type>;

    using omp_parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<omp_traits_type>;
    using seq_parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<seq_traits_type>;
    using omp_parameter_type      = typename omp_parameter_list_type::parameter_type;
    using seq_parameter_type      = typename seq_parameter_list_type::parameter_type;

    using omp_partition_type      = mjolnir::VerletList<omp_traits_type, potential_type>;
    using seq_partition_type      = mjolnir::VerletList<seq_traits_type, potential_type>;

    // those actually does not depends on traits.
    using base_kind            = typename omp_parameter_list_type::base_kind;
    using ignore_group_type    = typename omp_parameter_list_type::ignore_group_type;
    using ignore_molecule_type = typename omp_parameter_list_type::ignore_molecule_type;

    constexpr real_type pi = mjolnir::math::constants<real_type>::pi();

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal() / 1000.0) *
                          unit_type::avogadro_constant());
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal()) /
                            unit_type::avogadro_constant());
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom());

        phys_type::set_m_to_length(unit_type::m_to_angstrom());
        phys_type::set_length_to_m(unit_type::angstrom_to_m());

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom(), 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m(), 3));

        phys_type::set_length_unit(length);
    }

    //       Si   Bi   Bj   Sj
    //  5'    o -- o===o -- o     3'
    //  ^    /      \ /      \    |
    //  | P o        X        o P |
    //  |    \      / \      /    v
    //  3'    o -- o===o -- o     5'
    //       Bj_next   Bj_next
    //
    // Required test cases are the combination of the following conditions
    //
    // CrossStacking
    // theta3:
    //  1. theta < pi/2K
    //  2. pi/2K < theta < pi/K
    //  3. pi/K  < theta
    // thetaCS_i:
    //  1. theta < pi/2K
    //  2. pi/2K < theta < pi/K
    //  3. pi/K  < theta
    // thetaCS_j:
    //  1. theta < pi/2K
    //  2. pi/2K < theta < pi/K
    //  3. pi/K  < theta
    // thetaCS_j:
    //  1. theta < pi/2K
    //  2. pi/2K < theta < pi/K
    //  3. pi/K  < theta
    // ri_jnext
    //  1. (r < r0)
    //  2. (r0 < r)
    // rj_inext
    //  1. (r < r0)
    //  2. (r0 < r)
    //
    // Oh, okay, we have 3^4 * 2^2 = 324 test configurations.
    //
    // Also, this interaction is a sum of two, BasePairing and CrossStacking.
    // This makes numeric differentiation unstable.

    std::mt19937 rng(123456789);

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        const std::vector<base_kind> bases{base_kind::A, base_kind::T, base_kind::G, base_kind::G};
        //        0    1   9    8
        //  5'    o -- o===o -- o     3'
        //  ^    /      \ /      \    |
        //  | 2 o        X        o 7 |
        //  |    \      / \      /    v
        //  3'    o -- o===o -- o     5'
        //        3    4   6    5

        constexpr auto invalid = potential_type::invalid();

        omp_parameter_list_type omp_parameter_list(
            mjolnir::ThreeSPN2BaseBaseGlobalPotentialParameter<real_type>{}, {
                {1, omp_parameter_type{bases.at(0), 0, 0, 4, invalid}},
                {9, omp_parameter_type{bases.at(1), 3, 8, invalid, 6}},
                {4, omp_parameter_type{bases.at(2), 1, 3, invalid, 1}},
                {6, omp_parameter_type{bases.at(3), 2, 5, 9, invalid}},
            }, {}, ignore_molecule_type("Nothing"), ignore_group_type({})
        );
        seq_parameter_list_type seq_parameter_list(
            mjolnir::ThreeSPN2BaseBaseGlobalPotentialParameter<real_type>{}, {
                {1, seq_parameter_type{bases.at(0), 0, 0, 4, invalid}},
                {9, seq_parameter_type{bases.at(1), 3, 8, invalid, 6}},
                {4, seq_parameter_type{bases.at(2), 1, 3, invalid, 1}},
                {6, seq_parameter_type{bases.at(3), 2, 5, 9, invalid}},
            }, {}, ignore_molecule_type("Nothing"), ignore_group_type({})
        );

        omp_interaction_type omp_interaction(potential_type{},
            omp_parameter_list_type(omp_parameter_list),
            mjolnir::SpatialPartition<omp_traits_type, potential_type>(
                mjolnir::make_unique<omp_partition_type>()));

        seq_interaction_type seq_interaction(potential_type{},
            seq_parameter_list_type(seq_parameter_list),
            mjolnir::SpatialPartition<seq_traits_type, potential_type>(
                mjolnir::make_unique<seq_partition_type>()));

        omp_system_type omp_sys(10, boundary_type{});
        seq_system_type seq_sys(10, boundary_type{});

        test::clear_everything(omp_sys);
        test::clear_everything(seq_sys);
        topology_type topol(10);

        topol.add_connection(0, 1, "bond");
        topol.add_connection(0, 2, "bond");
        topol.add_connection(2, 3, "bond");
        topol.add_connection(3, 4, "bond");

        topol.add_connection(5, 6, "bond");
        topol.add_connection(5, 7, "bond");
        topol.add_connection(7, 8, "bond");
        topol.add_connection(8, 9, "bond");

        topol.add_connection(0, 2, "next_nucl");
        topol.add_connection(1, 2, "next_nucl");
        topol.add_connection(0, 3, "next_nucl");
        topol.add_connection(1, 3, "next_nucl");
        topol.add_connection(0, 4, "next_nucl");
        topol.add_connection(1, 4, "next_nucl");

        topol.add_connection(5, 7, "next_nucl");
        topol.add_connection(6, 7, "next_nucl");
        topol.add_connection(5, 8, "next_nucl");
        topol.add_connection(6, 8, "next_nucl");
        topol.add_connection(5, 9, "next_nucl");
        topol.add_connection(6, 9, "next_nucl");

        topol.construct_molecules();

        omp_interaction.initialize(omp_sys, topol);
        seq_interaction.initialize(seq_sys, topol);

        // paramter_lists in the interactions will be initialized via
        // interaction.initialize(), but parameter_list here will not.
        potential_type pot;
        omp_parameter_list.initialize(omp_sys, topol, pot);

        const auto bp_kind      = omp_parameter_list.bp_kind (bases.at(0), bases.at(1));
        const auto cs_i_kind    = omp_parameter_list.cs5_kind(bases.at(0), bases.at(3));
        const auto cs_j_kind    = omp_parameter_list.cs3_kind(bases.at(1), bases.at(2));

        const auto rbp          = omp_parameter_list.r0(bp_kind) + 0.2;
        const auto rcs_i_0      = omp_parameter_list.r0(cs_i_kind);
        const auto rcs_j_0      = omp_parameter_list.r0(cs_j_kind);
        const auto theta1       = omp_parameter_list.theta1_0(bp_kind);
        const auto theta2       = omp_parameter_list.theta2_0(bp_kind);
        const auto theta3_0     = omp_parameter_list.theta3_0(bp_kind);
        const auto thetaCS_i_0  = omp_parameter_list.thetaCS_0(cs_i_kind);
        const auto thetaCS_j_0  = omp_parameter_list.thetaCS_0(cs_j_kind);

        const auto pi_over_K_BP = omp_parameter_list.pi_over_K_BP();
        const auto pi_over_K_CS = omp_parameter_list.pi_over_K_CS();

        for(const auto rcsi : {rcs_i_0 - 0.2, rcs_i_0 + 0.5})
        {
        for(const auto rcsj : {rcs_j_0 - 0.2, rcs_j_0 + 0.5})
        {
        for(const auto theta3 : {theta3_0 + pi_over_K_BP * 0.1, theta3_0 + pi_over_K_BP * 0.6, theta3_0 + pi_over_K_BP * 1.1})
        {
        for(const auto thetaCS_i : {thetaCS_i_0 - pi_over_K_CS * 0.1, thetaCS_i_0 - pi_over_K_CS * 0.6, thetaCS_i_0 - pi_over_K_CS * 1.1})
        {
        for(const auto thetaCS_j : {thetaCS_j_0 + pi_over_K_CS * 0.1, thetaCS_j_0 + pi_over_K_CS * 0.6, thetaCS_j_0 + pi_over_K_CS * 1.1})
        {
            BOOST_TEST_REQUIRE(0.0 <= theta3   );BOOST_TEST_REQUIRE(theta3    <= pi);
            BOOST_TEST_REQUIRE(0.0 <= thetaCS_i);BOOST_TEST_REQUIRE(thetaCS_i <= pi);
            BOOST_TEST_REQUIRE(0.0 <= thetaCS_j);BOOST_TEST_REQUIRE(thetaCS_j <= pi);

            // generate positions...
            //        0    1   9    8
            //  5'    o -- o===o -- o     3'
            //  ^    /      \ /      \    |
            //  | 2 o        X        o 7 |
            //  |    \      / \      /    v
            //  3'    o -- o===o -- o     5'
            //        3    4   6    5

            omp_sys.position(1) = coordinate_type(0.0, 0.0, 0.0);
            omp_sys.position(9) = coordinate_type(0.0, 0.0, 0.0);

            omp_sys.position(0) = coordinate_type(4.0 * std::cos(theta1),
                                                  4.0 * std::sin(theta1), 0.0);
            omp_sys.position(8) = coordinate_type(4.0 * std::cos(pi - theta2),
                                                  4.0 * std::sin(pi - theta2), 0.0);
            // rotate around x axis
            // to make angle between 0->1 and 8->9 theta3
            {
                real_type best_phi = 0.0;
                real_type dtheta   = std::numeric_limits<real_type>::max();
                const auto bck = omp_sys;
                for(int i=-999; i<1000; ++i)
                {
                    const real_type phi = pi * 0.001 * i;
                    const matrix33_type rot(
                       1.0,           0.0,           0.0,
                       0.0, std::cos(phi), -std::sin(phi),
                       0.0, std::sin(phi),  std::cos(phi));

                    omp_sys.position(8) = rot * omp_sys.position(8);
                    omp_sys.position(9) = rot * omp_sys.position(9);

                    const auto v01 = omp_sys.position(1) - omp_sys.position(0);
                    const auto v89 = omp_sys.position(9) - omp_sys.position(8);
                    const auto theta_0189 = std::acos(mjolnir::math::clamp<real_type>(
                         mjolnir::math::dot_product(v01, v89) /
                        (mjolnir::math::length(v01) * mjolnir::math::length(v89)),
                        -1.0, 1.0));

                    BOOST_TEST_MESSAGE("dtheta = " << dtheta << ", current configuration = " << std::abs(theta_0189 - theta3));

                    if(std::abs(theta_0189 - theta3) < dtheta)
                    {
                        dtheta   = std::abs(theta_0189 - theta3);
                        best_phi = phi;
                    }
                    omp_sys = bck;
                }

                BOOST_TEST_MESSAGE("best_phi = " << best_phi << ", dtheta = " << dtheta);

                {
                    const matrix33_type rot(
                            1.0,                0.0,                 0.0,
                            0.0, std::cos(best_phi), -std::sin(best_phi),
                            0.0, std::sin(best_phi),  std::cos(best_phi));
                    omp_sys.position(8) = rot * omp_sys.position(8);
                    omp_sys.position(9) = rot * omp_sys.position(9);
                }

                mjolnir::math::X(omp_sys.position(4)) += rbp;
                mjolnir::math::X(omp_sys.position(8)) += rbp;
                mjolnir::math::X(omp_sys.position(9)) += rbp;
            }

            {
                const auto v01 = omp_sys.position(1) - omp_sys.position(0);
                const auto v89 = omp_sys.position(9) - omp_sys.position(8);
                const auto n4 = mjolnir::math::cross_product(coordinate_type(0.0, 0.0, 1.0), v89);
                const auto n6 = mjolnir::math::cross_product(coordinate_type(0.0, 0.0, 1.0), v01);

                const auto n4_reg = n4 / mjolnir::math::length(n4);
                const auto n6_reg = n6 / mjolnir::math::length(n6);

                coordinate_type v16 = v01 / mjolnir::math::length(v01);
                coordinate_type v94 = v89 / mjolnir::math::length(v89);
                {
                    const auto cos_theta = std::cos(pi - thetaCS_i);
                    const auto sin_theta = std::sin(pi - thetaCS_i);

                    BOOST_TEST_REQUIRE(mjolnir::math::length(n6_reg) == 1.0, boost::test_tools::tolerance(1e-8));
                    const auto nx = mjolnir::math::X(n6_reg);
                    const auto ny = mjolnir::math::Y(n6_reg);
                    const auto nz = mjolnir::math::Z(n6_reg);
                    const matrix33_type rot_16(
                        cos_theta + nx*nx*(1-cos_theta),    nx*ny*(1-cos_theta) - nz*sin_theta, nx*nz*(1-cos_theta) + ny*sin_theta,
                        ny*nx*(1-cos_theta) + nz*sin_theta, cos_theta + ny*ny*(1-cos_theta),    ny*nz*(1-cos_theta) - nx*sin_theta,
                        nz*nx*(1-cos_theta) - ny*sin_theta, nz*ny*(1-cos_theta) + nx*sin_theta, cos_theta + nz*nz*(1-cos_theta));
                    BOOST_TEST_REQUIRE(mjolnir::math::length(v16) == 1.0, boost::test_tools::tolerance(1e-8));
                    v16 = rot_16 * v16;
                    BOOST_TEST_REQUIRE(mjolnir::math::length(v16) == 1.0, boost::test_tools::tolerance(1e-8));
                }
                {
                    const auto cos_theta = std::cos(pi - thetaCS_j);
                    const auto sin_theta = std::sin(pi - thetaCS_j);

                    BOOST_TEST_REQUIRE(mjolnir::math::length(n4_reg) == 1.0, boost::test_tools::tolerance(1e-8));
                    const auto nx = mjolnir::math::X(n4_reg);
                    const auto ny = mjolnir::math::Y(n4_reg);
                    const auto nz = mjolnir::math::Z(n4_reg);

                    const matrix33_type rot_94(
                        cos_theta + nx*nx*(1-cos_theta),    nx*ny*(1-cos_theta) - nz*sin_theta, nx*nz*(1-cos_theta) + ny*sin_theta,
                        ny*nx*(1-cos_theta) + nz*sin_theta, cos_theta + ny*ny*(1-cos_theta),    ny*nz*(1-cos_theta) - nx*sin_theta,
                        nz*nx*(1-cos_theta) - ny*sin_theta, nz*ny*(1-cos_theta) + nx*sin_theta, cos_theta + nz*nz*(1-cos_theta));

                    BOOST_TEST_REQUIRE(mjolnir::math::length(v94) == 1.0, boost::test_tools::tolerance(1e-8));
                    v94 = rot_94 * v94;
                    BOOST_TEST_REQUIRE(mjolnir::math::length(v94) == 1.0, boost::test_tools::tolerance(1e-8));
                }
                omp_sys.position(4) = omp_sys.position(9) + v94 * rcsj;
                omp_sys.position(6) = omp_sys.position(1) + v16 * rcsi;
            }

            // check the configurations are okay
            {
                const auto v19 = omp_sys.position(9) - omp_sys.position(1);
                BOOST_TEST_REQUIRE(mjolnir::math::length(v19) == rbp,
                           boost::test_tools::tolerance(0.01));

                const auto v01 = omp_sys.position(1) - omp_sys.position(0);
                const auto v89 = omp_sys.position(9) - omp_sys.position(8);
                {
                    const auto dot = mjolnir::math::dot_product(v01, v89);
                    const auto cos_theta = dot /
                        (mjolnir::math::length(v01) * mjolnir::math::length(v89));
                    const auto theta = std::acos(cos_theta);
                    BOOST_TEST_REQUIRE(theta == theta3, boost::test_tools::tolerance(0.01));
                }

                {
                    const auto dot = mjolnir::math::dot_product(-v01, v19);
                    const auto cos_theta = dot /
                        (mjolnir::math::length(v01) * mjolnir::math::length(v19));
                    const auto theta = std::acos(cos_theta);
                    BOOST_TEST_REQUIRE(theta == theta1, boost::test_tools::tolerance(0.01));
                }
                {
                    const auto dot = mjolnir::math::dot_product(v89, v19);
                    const auto cos_theta = dot /
                        (mjolnir::math::length(v89) * mjolnir::math::length(v19));
                    const auto theta = std::acos(cos_theta);
                    BOOST_TEST_REQUIRE(theta == theta2, boost::test_tools::tolerance(0.01));
                }

                const auto v61 = omp_sys.position(1) - omp_sys.position(6);
                BOOST_TEST_REQUIRE(mjolnir::math::length(v61) == rcsi,
                                   boost::test_tools::tolerance(0.01));
                {
                    const auto dot = mjolnir::math::dot_product(
                            v01 / mjolnir::math::length(v01),
                            v61 / mjolnir::math::length(v61));
                    const auto theta = std::acos(dot);
                    BOOST_TEST_REQUIRE(theta == thetaCS_i,
                            boost::test_tools::tolerance(0.01));
                }
                const auto v49 = omp_sys.position(9) - omp_sys.position(4);
                BOOST_TEST_REQUIRE(mjolnir::math::length(v49) == rcsj,
                           boost::test_tools::tolerance(0.01));
                {
                    const auto dot = mjolnir::math::dot_product(
                            v89 / mjolnir::math::length(v89),
                            v49 / mjolnir::math::length(v49));
                    const auto theta = std::acos(dot);
                    BOOST_TEST_REQUIRE(theta == thetaCS_j, boost::test_tools::tolerance(0.01));
                }
            }

            // positions generated!
            // ... and then rotate random direction to remove special axis

            test::apply_random_rotation(omp_sys, rng);
            test::apply_random_perturbation(omp_sys, rng, 0.01);
            test::clear_force(omp_sys);
            test::clear_force(seq_sys);

            for(std::size_t idx=0; idx<10; ++idx)
            {
                seq_sys.position(idx) = omp_sys.position(idx);
            }

            test::check_force_consistency              (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
            test::check_force_and_energy_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
            test::check_force_and_virial_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
            test::check_force_energy_virial_consistency(omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
            test::check_energy_consistency             (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);

        } // thetaCS_j
        } // thetaCS_i
        } // theta3
        } // rcsj
        } // rcsi
    }
}
