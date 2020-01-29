#define BOOST_TEST_MODULE "test_3spn2_base_base_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteraction.hpp>
#include <mjolnir/input/read_units.hpp>
#include <random>

using parameters_to_test = std::tuple<
    mjolnir::ThreeSPN2BaseBaseGlobalPotentialParameter<double>,
    mjolnir::ThreeSPN2CBaseBaseGlobalPotentialParameter<double>
>;

BOOST_AUTO_TEST_CASE_TEMPLATE(ThreeSPN2BasePairIntearction_numerical_diff,
        ParameterSet, parameters_to_test)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;

    using potential_type      = mjolnir::ThreeSPN2BaseBaseInteractionPotential<traits_type>;
    using base_kind           = typename potential_type::base_kind;
    using parameter_type      = typename potential_type::parameter_type;
    using partition_type      = mjolnir::VerletList<traits_type, potential_type>;

    using interaction_type     = mjolnir::ThreeSPN2BaseBaseInteraction<traits_type>;
    using ignore_group_type    = typename potential_type::ignore_group_type;
    using ignore_molecule_type = typename potential_type::ignore_molecule_type;

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
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    for(const auto bases : {
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
        constexpr auto invalid = potential_type::invalid();

        potential_type potential(ParameterSet{}, {
                {1, parameter_type{bases.at(0), 0, 0, invalid, invalid}},
                {3, parameter_type{bases.at(1), 1, 2, invalid, invalid}}
            }, {}, ignore_molecule_type("Nothing"), ignore_group_type({})
        );

        interaction_type interaction(potential_type(potential),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));

        system_type sys(4, boundary_type{});

        sys.topology().add_connection(0, 1, "bond");
        sys.topology().add_connection(2, 3, "bond");
        sys.topology().construct_molecules();

        for(std::size_t i=0; i<4; ++i)
        {
            sys.mass(i)  = 1.0;
            sys.rmass(i) = 1.0;
            sys.position(i) = coord_type(0.0, 0.0, 0.0);
            sys.velocity(i) = coord_type(0.0, 0.0, 0.0);
            sys.force(i)    = coord_type(0.0, 0.0, 0.0);
            sys.name(i)  = "DNA";
            sys.group(i) = "DNA";
        }
        potential  .initialize(sys, sys.topology());
        interaction.initialize(sys, sys.topology());

        const auto bp_kind      = potential.bp_kind (bases.at(0), bases.at(1));
        const auto rbp_0        = potential.r0(bp_kind);
        const auto theta1_0     = potential.theta1_0(bp_kind);
        const auto theta2_0     = potential.theta2_0(bp_kind);
        const auto pi_over_K_BP = potential.pi_over_K_BP();

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

            // rotate in random direction
            {
                using matrix33_type = typename traits_type::matrix33_type;

                const auto rot_x = uni(mt) * pi;
                const auto rot_y = uni(mt) * pi;
                const auto rot_z = uni(mt) * pi;

                matrix33_type rotm_x(1.0,             0.0,              0.0,
                                     0.0, std::cos(rot_x), -std::sin(rot_x),
                                     0.0, std::sin(rot_x),  std::cos(rot_x));
                matrix33_type rotm_y( std::cos(rot_y), 0.0,  std::sin(rot_y),
                                                  0.0, 1.0,              0.0,
                                     -std::sin(rot_y), 0.0,  std::cos(rot_y));
                matrix33_type rotm_z(std::cos(rot_z), -std::sin(rot_z), 0.0,
                                     std::sin(rot_z),  std::cos(rot_z), 0.0,
                                                 0.0,              0.0, 1.0);

                const matrix33_type rotm = rotm_x * rotm_y * rotm_z;
                for(std::size_t idx=0; idx<sys.size(); ++idx)
                {
                    sys.position(idx) = rotm * sys.position(idx);
                }
            }
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

            for(std::size_t idx=0; idx<sys.size(); ++idx)
            {
                sys.position(idx) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
                sys.force(idx)     = coord_type(0.0, 0.0, 0.0);
            }
            const system_type init = sys;

            constexpr real_type tol = 1e-4;
            constexpr real_type dr  = 1e-5;
            for(std::size_t idx=0; idx<sys.size(); ++idx)
            {
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::X(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::X(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::X(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::Y(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::Y(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::Y(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::Z(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::Z(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::Z(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
            }
        } // theta2
        } // theta1
        } // rbp
    }

}

BOOST_AUTO_TEST_CASE_TEMPLATE(ThreeSPN2CrossStackingIntearction_numerical_diff,
        ParameterSet, parameters_to_test)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;

    using potential_type    = mjolnir::ThreeSPN2BaseBaseInteractionPotential<traits_type>;
    using base_kind         = typename potential_type::base_kind;
    using parameter_type    = typename potential_type::parameter_type;
    using partition_type      = mjolnir::VerletList<traits_type, potential_type>;

    using interaction_type     = mjolnir::ThreeSPN2BaseBaseInteraction<traits_type>;
    using ignore_group_type    = typename potential_type::ignore_group_type;
    using ignore_molecule_type = typename potential_type::ignore_molecule_type;

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
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    for(const auto bases : {
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

        potential_type potential(ParameterSet{}, {
                {1, parameter_type{bases.at(0), 0, 0, 4, invalid}},
                {9, parameter_type{bases.at(1), 3, 8, invalid, 6}},
                {4, parameter_type{bases.at(2), 1, 3, invalid, 1}},
                {6, parameter_type{bases.at(3), 2, 5, 9, invalid}},
            }, {}, ignore_molecule_type("Nothing"), ignore_group_type({})
        );

        interaction_type interaction(potential_type(potential),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));

        system_type sys(10, boundary_type{});

        sys.topology().add_connection(0, 1, "bond");
        sys.topology().add_connection(0, 2, "bond");
        sys.topology().add_connection(2, 3, "bond");
        sys.topology().add_connection(3, 4, "bond");

        sys.topology().add_connection(5, 6, "bond");
        sys.topology().add_connection(5, 7, "bond");
        sys.topology().add_connection(7, 8, "bond");
        sys.topology().add_connection(8, 9, "bond");

        sys.topology().add_connection(0, 2, "next_nucl");
        sys.topology().add_connection(1, 2, "next_nucl");
        sys.topology().add_connection(0, 3, "next_nucl");
        sys.topology().add_connection(1, 3, "next_nucl");
        sys.topology().add_connection(0, 4, "next_nucl");
        sys.topology().add_connection(1, 4, "next_nucl");

        sys.topology().add_connection(5, 7, "next_nucl");
        sys.topology().add_connection(6, 7, "next_nucl");
        sys.topology().add_connection(5, 8, "next_nucl");
        sys.topology().add_connection(6, 8, "next_nucl");
        sys.topology().add_connection(5, 9, "next_nucl");
        sys.topology().add_connection(6, 9, "next_nucl");

        sys.topology().construct_molecules();

        for(std::size_t i=0; i<10; ++i)
        {
            sys.mass(i)  = 1.0;
            sys.rmass(i) = 1.0;
            sys.position(i) = coord_type(0.0, 0.0, 0.0);
            sys.velocity(i) = coord_type(0.0, 0.0, 0.0);
            sys.force(i)    = coord_type(0.0, 0.0, 0.0);
            sys.name(i)  = "DNA";
            sys.group(i) = "DNA";
        }
        potential  .initialize(sys, sys.topology());
        interaction.initialize(sys, sys.topology());

        const auto bp_kind      = potential.bp_kind (bases.at(0), bases.at(1));
        const auto cs_i_kind    = potential.cs5_kind(bases.at(0), bases.at(3));
        const auto cs_j_kind    = potential.cs3_kind(bases.at(1), bases.at(2));

        const auto rbp          = potential.r0(bp_kind) + 0.2;
        const auto rcs_i_0      = potential.r0(cs_i_kind);
        const auto rcs_j_0      = potential.r0(cs_j_kind);
        const auto theta1       = potential.theta1_0(bp_kind);
        const auto theta2       = potential.theta2_0(bp_kind);
        const auto theta3_0     = potential.theta3_0(bp_kind);
        const auto thetaCS_i_0  = potential.thetaCS_0(cs_i_kind);
        const auto thetaCS_j_0  = potential.thetaCS_0(cs_j_kind);

        const auto pi_over_K_BP = potential.pi_over_K_BP();
        const auto pi_over_K_CS = potential.pi_over_K_CS();

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
            {
                using matrix33_type = typename traits_type::matrix33_type;

                const auto rot_x = uni(mt) * pi;
                const auto rot_y = uni(mt) * pi;
                const auto rot_z = uni(mt) * pi;

                matrix33_type rotm_x(1.0,             0.0,              0.0,
                                     0.0, std::cos(rot_x), -std::sin(rot_x),
                                     0.0, std::sin(rot_x),  std::cos(rot_x));
                matrix33_type rotm_y( std::cos(rot_y), 0.0,  std::sin(rot_y),
                                                  0.0, 1.0,              0.0,
                                     -std::sin(rot_y), 0.0,  std::cos(rot_y));
                matrix33_type rotm_z(std::cos(rot_z), -std::sin(rot_z), 0.0,
                                     std::sin(rot_z),  std::cos(rot_z), 0.0,
                                                 0.0,              0.0, 1.0);

                const matrix33_type rotm = rotm_x * rotm_y * rotm_z;
                for(std::size_t idx=0; idx<sys.size(); ++idx)
                {
                    sys.position(idx) = rotm * sys.position(idx);
                }
            }

            // add perturbation and reset force
            for(std::size_t idx=0; idx<10; ++idx)
            {
                sys.position(idx) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
                sys.force(idx)     = coord_type(0.0, 0.0, 0.0);
            }
            const system_type init = sys;

            constexpr real_type tol = 1e-4;
            constexpr real_type dr  = 1e-5;
            for(const std::size_t idx : {0, 1, 8, 9, 4, 6})
            {
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::X(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::X(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::X(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::Y(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::Y(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::Y(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::Z(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::Z(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::Z(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
            }
        } // thetaCS_j
        } // thetaCS_i
        } // theta3
        } // rcsj
        } // rcsi
    } // bp_kind
}
