#define BOOST_TEST_MODULE "test_3spn2_base_stacking_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/ThreeSPN2BaseStackingInteraction.hpp>
#include <mjolnir/input/read_units.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(ThreeSPN2BaseStackingInteraction_numerical_diff)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_stacking_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using interaction_type = mjolnir::ThreeSPN2BaseStackingInteraction<traits_type>;
    using potential_type   = mjolnir::ThreeSPN2BaseStackingPotential<real_type>;
    using base_stack_kind  = typename potential_type::base_stack_kind;

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

    //        SBi
    //     Si --> Bi
    //    /     `-^
    //   Pj theta | rij
    //    \       |
    //     Sj --- Bj
    //
    //  rij:
    //  1. (r < r0)
    //  2. (r0 < r)
    //  theta:
    //  1. theta < pi/2K
    //  2. pi/2K < theta < pi/K
    //  3. pi/K  < theta

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    constexpr real_type pi = mjolnir::math::constants<real_type>::pi;

    for(const auto bs_kind : {base_stack_kind::AA, base_stack_kind::AT,
                              base_stack_kind::AG, base_stack_kind::AC,
                              base_stack_kind::TA, base_stack_kind::TT,
                              base_stack_kind::TG, base_stack_kind::TC,
                              base_stack_kind::GA, base_stack_kind::GT,
                              base_stack_kind::GG, base_stack_kind::GC,
                              base_stack_kind::CA, base_stack_kind::CT,
                              base_stack_kind::CG, base_stack_kind::CC})
    {
        potential_type   potential{};
        interaction_type interaction("none",
                std::vector<std::pair<std::array<std::size_t, 3>, base_stack_kind>>{
                    { {{0, 1, 2}}, bs_kind }
                }, potential_type{});

        system_type sys(3, boundary_type{});

        sys.mass(0) = 1.0;
        sys.mass(1) = 1.0;
        sys.mass(2) = 1.0;

        sys.rmass(0) = 1.0;
        sys.rmass(1) = 1.0;
        sys.rmass(2) = 1.0;

        sys.position(0) = coord_type(0.0, 0.0, 0.0);
        sys.position(1) = coord_type(0.0, 0.0, 0.0);
        sys.position(2) = coord_type(0.0, 0.0, 0.0);

        sys.velocity(0) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(1) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(2) = coord_type(0.0, 0.0, 0.0);

        sys.force(0)    = coord_type(0.0, 0.0, 0.0);
        sys.force(1)    = coord_type(0.0, 0.0, 0.0);
        sys.force(2)    = coord_type(0.0, 0.0, 0.0);

        sys.name(0)  = "Si";
        sys.name(1)  = "Bi";
        sys.name(2)  = "Bj";
        sys.group(0) = "DNA";
        sys.group(1) = "DNA";
        sys.group(2) = "DNA";

        potential.initialize(sys);
        interaction.initialize(sys);

        const auto theta0    = potential.theta_0(bs_kind);
        const auto pi_over_K = potential.pi_over_K_BS();
        const auto theta0_1  = theta0 + 0.2 * pi_over_K; //         dtheta < pi/2K
        const auto theta0_2  = theta0 + 0.7 * pi_over_K; // pi/2K < dtheta < pi/K
        const auto theta0_3  = theta0 + 1.2 * pi_over_K; // pi/K  < dtheta
        const auto r0_1      = potential.r0(bs_kind) - 0.2;
        const auto r0_2      = potential.r0(bs_kind) + 0.5;

        for(const auto r : {r0_1, r0_2})
        {
        for(const auto theta : {theta0_1, theta0_2, theta0_3})
        {
            BOOST_TEST_MESSAGE("======================================");
            BOOST_TEST_MESSAGE("r = " << r << ", theta = " << theta);

        for(std::size_t i=0; i<100; ++i)
        {
            // generate particle configuration in the following way
            //    y
            // Bj ^
            //  \ | theta0
            // r0\|-.
            // ---o-----o--> x
            //  Bi      Si
            //

            sys.position(0) = coord_type(4.0, 0.0, 0.0); // Si
            sys.position(1) = coord_type(0.0, 0.0, 0.0); // Bi
            sys.position(2) = coord_type(r * std::cos(theta), r * std::sin(theta), 0.0); // Bj
            sys.force(0) = coord_type(0.0, 0.0, 0.0);
            sys.force(1) = coord_type(0.0, 0.0, 0.0);
            sys.force(2) = coord_type(0.0, 0.0, 0.0);

            // ... and then rotate random direction to remove special axis
            //
            // Do this thousand times with different random numbers!
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
                sys.position(0) = rotm * sys.position(0);
                sys.position(1) = rotm * sys.position(1);
                sys.position(2) = rotm * sys.position(2);

                BOOST_TEST(mjolnir::math::length(sys.position(0)) == 4.0,
                           boost::test_tools::tolerance(1e-6));
                BOOST_TEST(mjolnir::math::length(sys.position(1)) == 0.0,
                           boost::test_tools::tolerance(1e-6));
                BOOST_TEST(mjolnir::math::length(sys.position(2)) == r,
                           boost::test_tools::tolerance(1e-6));
            }

            sys.position(0) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.position(1) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.position(2) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));

            const system_type init = sys;

            BOOST_TEST(mjolnir::math::X(init.force(0)) == 0.0);
            BOOST_TEST(mjolnir::math::Y(init.force(0)) == 0.0);
            BOOST_TEST(mjolnir::math::Z(init.force(0)) == 0.0);

            BOOST_TEST(mjolnir::math::X(init.force(1)) == 0.0);
            BOOST_TEST(mjolnir::math::Y(init.force(1)) == 0.0);
            BOOST_TEST(mjolnir::math::Z(init.force(1)) == 0.0);

            BOOST_TEST(mjolnir::math::X(init.force(2)) == 0.0);
            BOOST_TEST(mjolnir::math::Y(init.force(2)) == 0.0);
            BOOST_TEST(mjolnir::math::Z(init.force(2)) == 0.0);

            constexpr real_type tol = 1e-4;
            constexpr real_type dr  = 1e-4;
            for(std::size_t idx=0; idx<3; ++idx)
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

                    BOOST_TEST(mjolnir::math::X(init.force(0)) == 0.0);
                    BOOST_TEST(mjolnir::math::Y(init.force(0)) == 0.0);
                    BOOST_TEST(mjolnir::math::Z(init.force(0)) == 0.0);

                    BOOST_TEST(mjolnir::math::X(init.force(1)) == 0.0);
                    BOOST_TEST(mjolnir::math::Y(init.force(1)) == 0.0);
                    BOOST_TEST(mjolnir::math::Z(init.force(1)) == 0.0);

                    BOOST_TEST(mjolnir::math::X(init.force(2)) == 0.0);
                    BOOST_TEST(mjolnir::math::Y(init.force(2)) == 0.0);
                    BOOST_TEST(mjolnir::math::Z(init.force(2)) == 0.0);

                    BOOST_TEST(mjolnir::math::X(sys.force(0)) == 0.0);
                    BOOST_TEST(mjolnir::math::Y(sys.force(0)) == 0.0);
                    BOOST_TEST(mjolnir::math::Z(sys.force(0)) == 0.0);

                    BOOST_TEST(mjolnir::math::X(sys.force(1)) == 0.0);
                    BOOST_TEST(mjolnir::math::Y(sys.force(1)) == 0.0);
                    BOOST_TEST(mjolnir::math::Z(sys.force(1)) == 0.0);

                    BOOST_TEST(mjolnir::math::X(sys.force(2)) == 0.0);
                    BOOST_TEST(mjolnir::math::Y(sys.force(2)) == 0.0);
                    BOOST_TEST(mjolnir::math::Z(sys.force(2)) == 0.0);

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
        } // theta
        } // r
    }
}
