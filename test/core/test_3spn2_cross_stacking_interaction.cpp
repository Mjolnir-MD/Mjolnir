#define BOOST_TEST_MODULE "test_3spn2_cross_stacking_interaction"

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
#include <mjolnir/forcefield/3SPN2/ThreeSPN2CrossStackingInteraction.hpp>
#include <mjolnir/input/read_units.hpp>
#include <random>

using parameters_to_test = std::tuple<
    mjolnir::ThreeSPN2CrossStackingGlobalPotentialParameter<double>,
    mjolnir::ThreeSPN2CCrossStackingGlobalPotentialParameter<double>
>;
BOOST_AUTO_TEST_CASE_TEMPLATE(ThreeSPN2CrossStackingIntearction_numerical_diff,
        ParameterSet, parameters_to_test)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_cross_stacking_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;
    using topology_type     = mjolnir::Topology;

    using interaction_type     = mjolnir::ThreeSPN2CrossStackingInteraction<traits_type>;

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

    std::mt19937 mt(123456789);

    for(const auto& bases : {
            std::vector<base_kind>{base_kind::A, base_kind::T, base_kind::G, base_kind::G},
            std::vector<base_kind>{base_kind::A, base_kind::T, base_kind::C, base_kind::C},
            std::vector<base_kind>{base_kind::T, base_kind::A, base_kind::G, base_kind::G},
            std::vector<base_kind>{base_kind::T, base_kind::A, base_kind::C, base_kind::C},
            std::vector<base_kind>{base_kind::C, base_kind::G, base_kind::A, base_kind::A},
            std::vector<base_kind>{base_kind::C, base_kind::G, base_kind::T, base_kind::T},
            std::vector<base_kind>{base_kind::G, base_kind::C, base_kind::A, base_kind::A},
            std::vector<base_kind>{base_kind::G, base_kind::C, base_kind::T, base_kind::T}
        })
    {
        //        0    1   9    8
        //  5'    o -- o===o -- o     3'
        //  ^    /      \ /      \    |
        //  | 2 o        X        o 7 |
        //  |    \      / \      /    v
        //  3'    o -- o===o -- o     5'
        //        3    4   6    5

        constexpr auto invalid = potential_type::invalid();

        potential_type potential;
        parameter_list_type parameter_list(ParameterSet{}, {
                {1, parameter_type{bases.at(0), 0, 4, invalid}},
                {9, parameter_type{bases.at(1), 8, invalid, 6}},
                {4, parameter_type{bases.at(2), 3, invalid, 1}},
                {6, parameter_type{bases.at(3), 5, 9, invalid}}
            }, {} /*{{"nucleotide", 3}}*/, ignore_molecule_type("Nothing"), ignore_group_type({})
        );
        mjolnir::ThreeSPN2BasePairParameterList<traits_type> bp_parameter_list(
            mjolnir::ThreeSPN2BasePairGlobalPotentialParameter<real_type>{}, {}, {},
            ignore_molecule_type("Nothing"), ignore_group_type({}));

        interaction_type interaction(potential_type{},
            parameter_list_type(parameter_list),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));

        system_type sys(10, boundary_type{});
        test::clear_everything(sys);

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

        parameter_list.initialize(sys, topol, potential);
        interaction.initialize(sys, topol);

        const auto bp_kind      = parameter_list.bp_kind(bases.at(0), bases.at(1));
        const auto cs_i_kind    = parameter_list.cs5_kind(bases.at(0), bases.at(3));
        const auto cs_j_kind    = parameter_list.cs3_kind(bases.at(1), bases.at(2));

        const auto rbp          = 6.0;
        const auto rcs_i_0      = parameter_list.r0(cs_i_kind);
        const auto rcs_j_0      = parameter_list.r0(cs_j_kind);
        const auto theta1       = bp_parameter_list.theta1_0(bp_kind);
        const auto theta2       = bp_parameter_list.theta2_0(bp_kind);
        const auto theta3_0     = parameter_list.theta3_0(bp_kind);
        const auto thetaCS_i_0  = parameter_list.thetaCS_0(cs_i_kind);
        const auto thetaCS_j_0  = parameter_list.thetaCS_0(cs_j_kind);

        const auto pi_over_K_BP = parameter_list.pi_over_K_BP();
        const auto pi_over_K_CS = parameter_list.pi_over_K_CS();

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
            BOOST_TEST_REQUIRE(0.0 <= theta3   ); BOOST_TEST_REQUIRE(theta3    <= pi);
            BOOST_TEST_REQUIRE(0.0 <= thetaCS_i); BOOST_TEST_REQUIRE(thetaCS_i <= pi);
            BOOST_TEST_REQUIRE(0.0 <= thetaCS_j); BOOST_TEST_REQUIRE(thetaCS_j <= pi);

            // generate positions...
            //        0    1   9    8
            //  5'    o -- o===o -- o     3'
            //  ^    /      \ /      \    |
            //  | 2 o        X        o 7 |
            //  |    \      / \      /    v
            //  3'    o -- o===o -- o     5'
            //        3    4   6    5

            sys.position(1) = coord_type(0.0, 0.0, 0.0);
            sys.position(9) = coord_type(0.0, 0.0, 0.0);

            sys.position(0) = coord_type(4.0 * std::cos(theta1),
                                         4.0 * std::sin(theta1), 0.0);
            sys.position(8) = coord_type(4.0 * std::cos(pi - theta2),
                                         4.0 * std::sin(pi - theta2), 0.0);
            // rotate around x axis
            // to make angle between 0->1 and 8->9 theta3
            {
                using matrix33_type = typename traits_type::matrix33_type;
                real_type best_phi = 0.0;
                real_type dtheta   = std::numeric_limits<real_type>::max();
                const auto bck = sys;
                for(int i=-999; i<1000; ++i)
                {
                    const real_type phi = pi * 0.001 * i;
                    const matrix33_type rot(
                       1.0,           0.0,           0.0,
                       0.0, std::cos(phi), -std::sin(phi),
                       0.0, std::sin(phi),  std::cos(phi));

                    sys.position(8) = rot * sys.position(8);
                    sys.position(9) = rot * sys.position(9);

                    const auto v01 = sys.position(1) - sys.position(0);
                    const auto v89 = sys.position(9) - sys.position(8);
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
                    sys = bck;
                }

                BOOST_TEST_MESSAGE("best_phi = " << best_phi << ", dtheta = " << dtheta);

                {
                    const matrix33_type rot(
                            1.0,                0.0,                 0.0,
                            0.0, std::cos(best_phi), -std::sin(best_phi),
                            0.0, std::sin(best_phi),  std::cos(best_phi));
                    sys.position(8) = rot * sys.position(8);
                    sys.position(9) = rot * sys.position(9);
                }

                mjolnir::math::X(sys.position(4)) += rbp;
                mjolnir::math::X(sys.position(8)) += rbp;
                mjolnir::math::X(sys.position(9)) += rbp;
            }

            {
                using matrix33_type = typename traits_type::matrix33_type;
                const auto v01 = sys.position(1) - sys.position(0);
                const auto v89 = sys.position(9) - sys.position(8);
                const auto n4 = mjolnir::math::cross_product(coord_type(0.0, 0.0, 1.0), v89);
                const auto n6 = mjolnir::math::cross_product(coord_type(0.0, 0.0, 1.0), v01);

                const auto n4_reg = n4 / mjolnir::math::length(n4);
                const auto n6_reg = n6 / mjolnir::math::length(n6);

                coord_type v16 = v01 / mjolnir::math::length(v01);
                coord_type v94 = v89 / mjolnir::math::length(v89);
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
                sys.position(4) = sys.position(9) + v94 * rcsj;
                sys.position(6) = sys.position(1) + v16 * rcsi;
            }

            // check the configurations are okay
            {
                const auto v19 = sys.position(9) - sys.position(1);
                BOOST_TEST_REQUIRE(mjolnir::math::length(v19) == rbp,
                           boost::test_tools::tolerance(0.01));

                const auto v01 = sys.position(1) - sys.position(0);
                const auto v89 = sys.position(9) - sys.position(8);
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

                const auto v61 = sys.position(1) - sys.position(6);
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
                const auto v49 = sys.position(9) - sys.position(4);
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

            test::apply_random_rotation(sys, mt);
            test::apply_random_perturbation(sys, mt, 0.01);
            test::clear_force(sys);

            constexpr real_type tol = 1e-4;
            constexpr real_type dr  = 1e-5;
            test::check_force (sys, interaction, tol, dr);

            // -----------------------------------------------------------------
            // check virial

            test::check_virial(sys, interaction, tol);
            test::check_force_and_virial(sys, interaction, tol);
            test::check_force_and_energy(sys, interaction, tol);
            test::check_force_energy_virial(sys, interaction, tol);

        } // thetaCS_j
        } // thetaCS_i
        } // theta3
        } // rcsj
        } // rcsi
    } // bp_kind
}
