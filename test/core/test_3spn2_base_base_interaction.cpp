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
#include <mjolnir/interaction/global/ThreeSPN2BaseBaseInteraction.hpp>
#include <mjolnir/input/read_units.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(ThreeSPN2BasePairIntearction_numerical_diff)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;

    using potential_type      = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;
    using base_kind           = typename potential_type::base_kind;
    using parameter_type      = typename potential_type::parameter_type;
    using pair_parameter_type = typename potential_type::pair_parameter_type;
    using partition_type      = mjolnir::VerletList<traits_type, pair_parameter_type>;

    using interaction_type  = mjolnir::ThreeSPN2BaseBaseInteraction<traits_type, partition_type>;
    using ignore_group_type = typename potential_type::ignore_group_type;

    constexpr real_type pi = mjolnir::math::constants<real_type>::pi;

    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal / 1000.0) *
                          unit_type::avogadro_constant);
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal) /
                            unit_type::avogadro_constant);
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom);

        phys_type::set_m_to_length(unit_type::m_to_angstrom);
        phys_type::set_length_to_m(unit_type::angstrom_to_m);

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom, 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m, 3));

        phys_type::set_length_unit("angstrom");
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

        ignore_group_type grp({});
        potential_type potential({
                {1, parameter_type{bases.at(0), 0, 0, invalid, invalid}},
                {3, parameter_type{bases.at(1), 1, 2, invalid, invalid}}
            }, std::move(grp)
        );

        interaction_type interaction(
                potential_type(potential), partition_type{});

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
        potential  .initialize(sys);
        interaction.initialize(sys);

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
        for(std::size_t i=0; i<100; ++i) // perturbation
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
            for(std::size_t idx=0; idx<sys.size(); ++idx)
            {
                sys.position(idx) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
                sys.force(idx)     = coord_type(0.0, 0.0, 0.0);
            }
            const system_type init = sys;

            constexpr real_type tol = 1e-3;
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

                    BOOST_TEST(-dE == dr * mjolnir::math::X(sys.force(idx)),
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

                    BOOST_TEST(-dE == dr * mjolnir::math::Y(sys.force(idx)),
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

                    BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
            }
        } // perturbation
        } // theta2
        } // theta1
        } // rbp
    }

}

BOOST_AUTO_TEST_CASE(ThreeSPN2CrossStackingIntearction_numerical_diff)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;

    using potential_type    = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;
    using base_kind         = typename potential_type::base_kind;
    using parameter_type    = typename potential_type::parameter_type;
    using pair_parameter_type = typename potential_type::pair_parameter_type;
    using partition_type      = mjolnir::VerletList<traits_type, pair_parameter_type>;

    using interaction_type  = mjolnir::ThreeSPN2BaseBaseInteraction<traits_type, partition_type>;
    using ignore_group_type = typename potential_type::ignore_group_type;

    constexpr real_type pi = mjolnir::math::constants<real_type>::pi;

    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal / 1000.0) *
                          unit_type::avogadro_constant);
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal) /
                            unit_type::avogadro_constant);
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom);

        phys_type::set_m_to_length(unit_type::m_to_angstrom);
        phys_type::set_length_to_m(unit_type::angstrom_to_m);

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom, 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m, 3));

        phys_type::set_length_unit("angstrom");
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

        potential_type potential({
                {1, parameter_type{bases.at(0), 0, 0, 4, invalid}},
                {9, parameter_type{bases.at(1), 3, 8, invalid, 6}},
                {4, parameter_type{bases.at(2), 1, 3, invalid, 1}},
                {6, parameter_type{bases.at(3), 2, 5, 9, invalid}},
            }, ignore_group_type({})
        );

        interaction_type interaction(
                potential_type(potential), partition_type{});

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
        potential  .initialize(sys);
        interaction.initialize(sys);

        const auto bp_kind      = potential.bp_kind (bases.at(0), bases.at(1));
        const auto cs_i_kind    = potential.cs5_kind(bases.at(0), bases.at(3));
        const auto cs_j_kind    = potential.cs3_kind(bases.at(1), bases.at(2));

        const auto rbp          = potential.r0(bp_kind) + 0.2;
        const auto rcs_i_0      = potential.r0(cs_i_kind);
        const auto rcs_j_0      = potential.r0(cs_j_kind);
        const auto theta3_0     = potential.theta3_0(bp_kind);
        const auto thetaCS_i_0  = potential.thetaCS_0(cs_i_kind);
        const auto thetaCS_j_0  = potential.thetaCS_0(cs_j_kind);

        const auto pi_over_K_BP = potential.pi_over_K_BP();
        const auto pi_over_K_CS = potential.pi_over_K_CS();

        for(const auto rcsi : {rcs_i_0 - 0.2, rcs_i_0 + 0.5})
        {
        for(const auto rcsj : {rcs_j_0 - 0.2, rcs_j_0 + 0.5})
        {
        for(const auto theta3 : {theta3_0 + pi_over_K_BP * 0.2, theta3_0 + pi_over_K_BP * 0.7, theta3_0 + pi_over_K_BP * 1.2})
        {
        for(const auto thetaCS_i : {thetaCS_i_0 - pi_over_K_CS * 0.1, thetaCS_i_0 - pi_over_K_CS * 0.6, thetaCS_i_0 - pi_over_K_CS * 1.1})
        {
        for(const auto thetaCS_j : {thetaCS_j_0 + pi_over_K_CS * 0.1, thetaCS_j_0 + pi_over_K_CS * 0.6, thetaCS_j_0 + pi_over_K_CS * 1.1})
        {

        BOOST_TEST_MESSAGE("========================================");
        BOOST_TEST_MESSAGE("rcs_i  = " << (rcsi > rcs_i_0 ? "large" : "small"));
        BOOST_TEST_MESSAGE("rcs_j  = " << (rcsj > rcs_j_0 ? "large" : "small"));
        if(theta3 == theta3_0 + pi_over_K_BP * 0.2)
        {
            BOOST_TEST_MESSAGE("theta3 = full");
        }
        else if(theta3 == theta3_0 + pi_over_K_BP * 0.7)
        {
            BOOST_TEST_MESSAGE("theta3 = partial");
        }
        else
        {
            BOOST_TEST_MESSAGE("theta3 = none");
        }
        if(thetaCS_i == thetaCS_i_0 - pi_over_K_CS * 0.1)
        {
            BOOST_TEST_MESSAGE("thetaCS_i = full");
        }
        else if(thetaCS_i == thetaCS_i_0 - pi_over_K_CS * 0.6)
        {
            BOOST_TEST_MESSAGE("thetaCS_i = partial");
        }
        else
        {
            BOOST_TEST_MESSAGE("thetaCS_i = none");
        }

        if(thetaCS_j == thetaCS_j_0 + pi_over_K_CS * 0.1)
        {
            BOOST_TEST_MESSAGE("thetaCS_j = full");
        }
        else if(thetaCS_j == thetaCS_j_0 + pi_over_K_CS * 0.6)
        {
            BOOST_TEST_MESSAGE("thetaCS_j = partial");
        }
        else
        {
            BOOST_TEST_MESSAGE("thetaCS_j = none");
        }

        for(std::size_t i=0; i<10; ++i) // perturbation
        {
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

            sys.position(0) = coord_type(-14.0, 0.0, 0.0);
            sys.position(8) = coord_type( 14.0, 0.0, 0.0);

            sys.position(4) = coord_type(-rcsj * std::cos(pi - thetaCS_j),
                                          0.0,
                                         -rcsj * std::sin(pi - thetaCS_j));
            sys.position(6) = coord_type(-rcsi * std::cos(thetaCS_i),
                                          0.0,
                                         -rcsi * std::sin(thetaCS_i));

            // rotate both around z axis
            // to make angle between 0->1 and 8->9 theta3
            {
                using matrix33_type = typename traits_type::matrix33_type;
                const auto rot_angle_i = -0.5 * (pi - theta3);
                const matrix33_type rot_i(
                        std::cos(rot_angle_i), -std::sin(rot_angle_i), 0.0,
                        std::sin(rot_angle_i),  std::cos(rot_angle_i), 0.0,
                        0.0,                   0.0,                    1.0);

                sys.position(0) = rot_i * sys.position(0);
                sys.position(1) = rot_i * sys.position(1);
                sys.position(5) = rot_i * sys.position(5);
                sys.position(6) = rot_i * sys.position(6);

                const auto rot_angle_j = 0.5 * (pi - theta3);
                const matrix33_type rot_j(
                        std::cos(rot_angle_j), -std::sin(rot_angle_j), 0.0,
                        std::sin(rot_angle_j),  std::cos(rot_angle_j), 0.0,
                        0.0,                                      0.0, 1.0);

                sys.position(8) = rot_j * sys.position(8);
                sys.position(9) = rot_j * sys.position(9);
                sys.position(3) = rot_j * sys.position(3);
                sys.position(4) = rot_j * sys.position(4);

                mjolnir::math::X(sys.position(8)) += rbp;
                mjolnir::math::X(sys.position(9)) += rbp;
                mjolnir::math::X(sys.position(4)) += rbp;
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
            constexpr real_type dr  = 1e-4;
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

                    BOOST_TEST(-dE == dr * mjolnir::math::X(sys.force(idx)),
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

                    BOOST_TEST(-dE == dr * mjolnir::math::Y(sys.force(idx)),
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

                    BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
            }
        } // configuration
        } // thetaCS_j
        } // thetaCS_i
        } // theta3
        } // rcsj
        } // rcsi
    } // bp_kind
}
