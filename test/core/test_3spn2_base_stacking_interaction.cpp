#define BOOST_TEST_MODULE "test_dihedral_angle_interaction"

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

#include <random>

BOOST_AUTO_TEST_CASE(DihedralAngleInteraction_numerical_diff)
{
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using interaction_type = mjolnir::ThreeSPN2BaseStackingInteraction<traits_type>;
    using potential_type   = mjolnir::ThreeSPN2BaseStackingPotential<real_type>;
    using base_stack_kind  = typename potential_type::base_stack_kind;

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
        const potential_type   potential{};
        const interaction_type interaction("none",
                std::vector<std::pair<std::array<std::size_t, 3>, base_stack_kind>>{
                    { {{0, 1, 2}}, bs_kind }
                }, potential_type{});

        system_type sys(3, boundary_type{});

        sys.at(0).mass = 1.0;
        sys.at(1).mass = 1.0;
        sys.at(2).mass = 1.0;

        sys.at(0).rmass = 1.0;
        sys.at(1).rmass = 1.0;
        sys.at(2).rmass = 1.0;

        sys.at(0).position = coord_type(0.0, 0.0, 0.0);
        sys.at(1).position = coord_type(0.0, 0.0, 0.0);
        sys.at(2).position = coord_type(0.0, 0.0, 0.0);

        sys.at(0).velocity = coord_type(0.0, 0.0, 0.0);
        sys.at(1).velocity = coord_type(0.0, 0.0, 0.0);
        sys.at(2).velocity = coord_type(0.0, 0.0, 0.0);

        sys.at(0).force    = coord_type(0.0, 0.0, 0.0);
        sys.at(1).force    = coord_type(0.0, 0.0, 0.0);
        sys.at(2).force    = coord_type(0.0, 0.0, 0.0);

        sys.at(0).name  = "Si";
        sys.at(1).name  = "Bi";
        sys.at(2).name  = "Bj";
        sys.at(0).group = "DNA";
        sys.at(1).group = "DNA";
        sys.at(2).group = "DNA";

        for(const auto r0 : {potential.r0(bs_kind) - 0.5,  // 2 test cases
                             potential.r0(bs_kind) + 0.5}) // to check both
        {
            const auto theta0    = potential.theta_0(bs_kind);
            const auto pi_over_K = potential.pi_over_K_BS();
            const auto theta0_1  = theta0 + 0.2 * pi_over_K; //         dtheta < pi/2K
            const auto theta0_2  = theta0 + 0.7 * pi_over_K; // pi/2K < dtheta < pi/K
            const auto theta0_3  = theta0 + 1.2 * pi_over_K; // pi/K  < dtheta

            for(const auto theta : {theta0_1, theta0_2, theta0_3})
            {
                for(std::size_t i=0; i<1000; ++i)
                {
                    // generate particle configuration in the following way
                    //    y
                    // Bj ^
                    //  \ | theta0
                    // r0\|-.
                    // ---o-----o--> x
                    //  Bi      Si
                    //

                    sys.at(0).position = coord_type(4.0, 0.0, 0.0); // Si
                    sys.at(1).position = coord_type(0.0, 0.0, 0.0); // Bi
                    sys.at(2).position = coord_type(r0 * std::cos(theta), r0 * std::sin(theta), 0.0); // Bj

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
                        sys.at(0).position = rotm * sys.at(0).position;
                        sys.at(1).position = rotm * sys.at(1).position;
                        sys.at(2).position = rotm * sys.at(2).position;
                    }

                    const auto init = sys;

                    constexpr real_type tol = 1e-5;
                    constexpr real_type dr  = 1e-5;
                    for(std::size_t idx=0; idx<3; ++idx)
                    {
                        {
                            // ----------------------------------------------------------------
                            // reset positions
                            sys = init;

                            // calc U(x-dx)
                            const auto E0 = interaction.calc_energy(sys);

                            mjolnir::math::X(sys.position(0)) += dr;

                            // calc F(x)
                            interaction.calc_force(sys);

                            mjolnir::math::X(sys.position(0)) += dr;

                            // calc U(x+dx)
                            const auto E1 = interaction.calc_energy(sys);

                            // central difference
                            const auto dE = (E1 - E0) * 0.5;

                            BOOST_TEST(-dE == dr * mjolnir::math::X(sys.force(0)),
                                       boost::test_tools::tolerance(tol));
                        }
                        {
                            // ----------------------------------------------------------------
                            // reset positions
                            sys = init;

                            // calc U(x-dx)
                            const auto E0 = interaction.calc_energy(sys);

                            mjolnir::math::Y(sys.position(0)) += dr;

                            // calc F(x)
                            interaction.calc_force(sys);

                            mjolnir::math::Y(sys.position(0)) += dr;

                            // calc U(x+dx)
                            const auto E1 = interaction.calc_energy(sys);

                            // central difference
                            const auto dE = (E1 - E0) * 0.5;

                            BOOST_TEST(-dE == dr * mjolnir::math::Y(sys.force(0)),
                                       boost::test_tools::tolerance(tol));
                        }
                        {
                            // ----------------------------------------------------------------
                            // reset positions
                            sys = init;

                            // calc U(x-dx)
                            const auto E0 = interaction.calc_energy(sys);

                            mjolnir::math::Z(sys.position(0)) += dr;

                            // calc F(x)
                            interaction.calc_force(sys);

                            mjolnir::math::Z(sys.position(0)) += dr;

                            // calc U(x+dx)
                            const auto E1 = interaction.calc_energy(sys);

                            // central difference
                            const auto dE = (E1 - E0) * 0.5;

                            BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(0)),
                                       boost::test_tools::tolerance(tol));
                        }
                    }
                }
            }
        }
    }
}
