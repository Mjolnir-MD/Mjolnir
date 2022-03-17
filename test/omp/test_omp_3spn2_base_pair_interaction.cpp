#define BOOST_TEST_MODULE "test_omp_3spn2_base_pair_interaction"

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
#include <mjolnir/omp/ThreeSPN2BasePairInteraction.hpp>

BOOST_AUTO_TEST_CASE(omp_ThreeSPN2BasePairIntearction_basepair)
{
    namespace test = mjolnir::test;

    constexpr double tol = 1e-4;
    mjolnir::LoggerManager::set_default_logger(
            "test_omp_3spn2_base_pair_interaction.log");

    using seq_traits_type  = mjolnir::SimulatorTraits      <double, mjolnir::UnlimitedBoundary>;
    using omp_traits_type  = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = typename omp_traits_type::real_type;
    using coordinate_type  = typename omp_traits_type::coordinate_type;
    using boundary_type    = typename omp_traits_type::boundary_type;
    using topology_type    = mjolnir::Topology;

    using potential_type       = mjolnir::ThreeSPN2BasePairPotential<real_type>;

    using omp_system_type      = mjolnir::System<omp_traits_type>;
    using omp_interaction_type = mjolnir::ThreeSPN2BasePairInteraction<omp_traits_type>;

    using seq_system_type      = mjolnir::System<seq_traits_type>;
    using seq_interaction_type = mjolnir::ThreeSPN2BasePairInteraction<seq_traits_type>;

    using omp_parameter_list_type = mjolnir::ThreeSPN2BasePairParameterList<omp_traits_type>;
    using seq_parameter_list_type = mjolnir::ThreeSPN2BasePairParameterList<seq_traits_type>;
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

        omp_parameter_list_type omp_parameter_list(
            mjolnir::ThreeSPN2BasePairGlobalPotentialParameter<real_type>{}, {
                {1, omp_parameter_type{bases.at(0), 0}},
                {3, omp_parameter_type{bases.at(1), 2}}
            }, {}, ignore_molecule_type("Nothing"), ignore_group_type({})
        );
        seq_parameter_list_type seq_parameter_list(
            mjolnir::ThreeSPN2BasePairGlobalPotentialParameter<real_type>{}, {
                {1, seq_parameter_type{bases.at(0), 0}},
                {3, seq_parameter_type{bases.at(1), 2}}
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
