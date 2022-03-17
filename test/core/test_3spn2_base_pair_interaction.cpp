#define BOOST_TEST_MODULE "test_3spn2_base_pair_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BasePairInteraction.hpp>
#include <mjolnir/input/read_units.hpp>
#include <random>

using parameters_to_test = std::tuple<
    mjolnir::ThreeSPN2BasePairGlobalPotentialParameter<double>,
    mjolnir::ThreeSPN2CBasePairGlobalPotentialParameter<double>
>;

BOOST_AUTO_TEST_CASE_TEMPLATE(ThreeSPN2BasePairIntearction_numerical_diff,
        ParameterSet, parameters_to_test)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_pair_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;
    using topology_type     = mjolnir::Topology;

    using interaction_type     = mjolnir::ThreeSPN2BasePairInteraction<traits_type>;

    using potential_type       = typename interaction_type::potential_type;
    using parameter_list_type  = typename interaction_type::parameter_list_type;
    using parameter_type       = typename parameter_list_type::parameter_type;
    using base_kind            = typename parameter_list_type::base_kind;
    using ignore_group_type    = typename parameter_list_type::ignore_group_type;
    using ignore_molecule_type = typename parameter_list_type::ignore_molecule_type;

    using partition_type       = mjolnir::VerletList<traits_type, potential_type>;

    constexpr real_type pi = mjolnir::math::constants<real_type>::pi();
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

    std::mt19937 mt(123456789);

    for(const auto& bases : {
            std::vector<base_kind>{base_kind::A, base_kind::T},
            std::vector<base_kind>{base_kind::T, base_kind::A},
            std::vector<base_kind>{base_kind::C, base_kind::G},
            std::vector<base_kind>{base_kind::G, base_kind::C}
            })
    {
        //  theta1    theta2
        //       |    |
        //  0 o  v    v o 2
        //     \-.   .-/
        //      o === o
        //     1       3
        //

        potential_type      potential;
        parameter_list_type parameter_list(ParameterSet{}, {
                {1, parameter_type{bases.at(0), 0}},
                {3, parameter_type{bases.at(1), 2}}
            }, {}, ignore_molecule_type("Nothing"), ignore_group_type({})
        );

        interaction_type interaction(potential_type{},
            parameter_list_type(parameter_list),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));

        system_type sys(4, boundary_type{});
        test::clear_everything(sys);

        topology_type topol(4);

        topol.add_connection(0, 1, "bond");
        topol.add_connection(2, 3, "bond");
        topol.construct_molecules();

        parameter_list.initialize(sys, topol, potential);
        interaction.initialize(sys, topol);

        const auto bp_kind      = parameter_list.bp_kind (bases.at(0), bases.at(1));
        const auto rbp_0        = parameter_list.r0(bp_kind);
        const auto theta1_0     = parameter_list.theta1_0(bp_kind);
        const auto theta2_0     = parameter_list.theta2_0(bp_kind);
        const auto pi_over_K_BP = parameter_list.pi_over_K_BP();

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
            sys.position(1) = coord_type(0.0, 0.0, 0.0);
            sys.position(3) = coord_type(rbp, 0.0, 0.0);

            sys.position(0) = coord_type(std::cos(theta1),          std::sin(theta1),    0.0);
            sys.position(2) = coord_type(std::cos(pi-theta2) + rbp, std::sin(pi-theta2), 0.0);

            test::apply_random_rotation(sys, mt);

            // check configuration is okay
            {
                const auto v13 = sys.position(3) - sys.position(1);
                BOOST_TEST_REQUIRE(mjolnir::math::length(v13) == rbp,
                                   boost::test_tools::tolerance(1e-3));

                const auto v10 = sys.position(0) - sys.position(1);
                const auto cos_013 = mjolnir::math::dot_product(v13, v10) /
                    (mjolnir::math::length(v13) * mjolnir::math::length(v10));
                BOOST_TEST_REQUIRE(std::acos(cos_013) == theta1,
                                   boost::test_tools::tolerance(1e-3));

                const auto v23 = sys.position(3) - sys.position(2);
                const auto cos_132 = mjolnir::math::dot_product(v13, v23) /
                    (mjolnir::math::length(v13) * mjolnir::math::length(v23));
                BOOST_TEST_REQUIRE(std::acos(cos_132) == theta2,
                                   boost::test_tools::tolerance(1e-3));
            }

            test::apply_random_perturbation(sys, mt, 0.01);
            test::clear_force(sys);

            constexpr real_type tol = 1e-4;
            constexpr real_type dr  = 1e-5;

            test::check_force (sys, interaction, tol, dr);
            test::check_virial(sys, interaction, tol);
            test::check_force_and_virial(sys, interaction, tol);
            test::check_force_and_energy(sys, interaction, tol);
            test::check_force_energy_virial(sys, interaction, tol);

        } // theta2
        } // theta1
        } // rbp
    }
}
