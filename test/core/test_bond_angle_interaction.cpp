#define BOOST_TEST_MODULE "test_bond_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/BondAngleInteraction.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(BondAngleInteraction_force)
{
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = traits_type::real_type;
    using coord_type      = traits_type::coordinate_type;
    using boundary_type   = traits_type::boundary_type;
    using system_type     = mjolnir::System<traits_type>;
    using harmonic_type   = mjolnir::HarmonicPotential<real_type>;
    using bond_angle_type = mjolnir::BondAngleInteraction<traits_type, harmonic_type>;

    constexpr real_type tol = 1e-7;

    auto normalize = [](const coord_type& v){return v / mjolnir::math::length(v);};

    const real_type k(1e0);
    const real_type native(mjolnir::math::constants<real_type>::pi * 2.0 / 3.0); // 120 degree

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    const coord_type pos1(1., 0., 0.);
    const coord_type pos2(0., 0., 0.);
    system_type sys(3, boundary_type{});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(2).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;
    sys.at(2).rmass = 1.0;

    sys.at(0).position = pos1;
    sys.at(1).position = pos2;
    sys.at(2).position = coord_type(0,0,0);
    sys.at(0).velocity = coord_type(0,0,0);
    sys.at(1).velocity = coord_type(0,0,0);
    sys.at(2).velocity = coord_type(0,0,0);
    sys.at(0).force    = coord_type(0,0,0);
    sys.at(1).force    = coord_type(0,0,0);
    sys.at(2).force    = coord_type(0,0,0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(2).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";
    sys.at(2).group = "NONE";

    const int N = 1800;
    const real_type dtheta = mjolnir::math::constants<real_type>::pi  / N;
    for(int i = 1; i < N; ++i)
    {
        BOOST_TEST(mjolnir::math::length(sys[0].position - pos1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[1].position - pos2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[0].velocity)        == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[1].velocity)        == 0.0, boost::test_tools::tolerance(tol));
        sys[0].force = coord_type(0,0,0);
        sys[1].force = coord_type(0,0,0);
        sys[2].force = coord_type(0,0,0);

        const real_type theta = i * dtheta;
        const coord_type pos3(std::cos(theta), std::sin(theta), 0e0);
        sys[2].position = pos3;

        const real_type deriv = potential.derivative(theta);
        const real_type coef = std::abs(deriv);

        interaction.calc_force(sys);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV/dtheta.
        BOOST_TEST(mjolnir::math::length(sys[0].position - sys[1].position) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[2].position - sys[1].position) == 1.0, boost::test_tools::tolerance(tol));

        const real_type force_strength1 = mjolnir::math::length(sys[0].force);
        const real_type force_strength3 = mjolnir::math::length(sys[2].force);
        if(i == 1200) // most stable point
        {
            BOOST_TEST(coef == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == 0.0, boost::test_tools::tolerance(tol));
        }
        else
        {
            BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == force_strength3, boost::test_tools::tolerance(tol));
        }

        // force applied to center particle is equal to sum of others
        const coord_type sum = sys[0].force + sys[1].force + sys[2].force;
        BOOST_TEST(mjolnir::math::length(sum) == 0.0, boost::test_tools::tolerance(tol));

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_TEST(std::abs(force_strength1) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(std::abs(force_strength3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i < 1200) // narrow
        {
            // perpendicular to radius vector
            const real_type dot1 = mjolnir::math::dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = mjolnir::math::dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., -1., 0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[0].force) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta + mjolnir::math::constants<real_type>::pi / 2.),
                    sin(theta + mjolnir::math::constants<real_type>::pi / 2.),
                    0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[2].force) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i > 1200) // extensive
        {
            const real_type dot1 = mjolnir::math::dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = mjolnir::math::dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., 1., 0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[0].force) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta - mjolnir::math::constants<real_type>::pi / 2.),
                    sin(theta - mjolnir::math::constants<real_type>::pi / 2.),
                    0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[2].force) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }

        // perpendicular to z axis
        BOOST_TEST(sys[0].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[1].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[2].force[2] == 0.0, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(BondAngleInteraction_numerical_diff)
{
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = traits_type::real_type;
    using coord_type      = traits_type::coordinate_type;
    using boundary_type   = traits_type::boundary_type;
    using system_type     = mjolnir::System<traits_type>;
    using harmonic_type   = mjolnir::HarmonicPotential<real_type>;
    using bond_angle_type = mjolnir::BondAngleInteraction<traits_type, harmonic_type>;

    const real_type k(10.0);
    const real_type native(mjolnir::math::constants<real_type>::pi / 3.0);

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    const coord_type pos1;
    const coord_type pos2;
    system_type sys(3, boundary_type{});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(2).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;
    sys.at(2).rmass = 1.0;

    sys.at(0).position = coord_type(1.0, 0.0, 0.0);
    sys.at(1).position = coord_type(0.0, 0.0, 1.0);
    sys.at(2).position = coord_type(1.0, 1.0, 1.0);
    sys.at(0).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(1).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(2).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(0).force    = coord_type(0.0, 0.0, 0.0);
    sys.at(1).force    = coord_type(0.0, 0.0, 0.0);
    sys.at(2).force    = coord_type(0.0, 0.0, 0.0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(2).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";
    sys.at(2).group = "NONE";

    const auto init = sys;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    constexpr real_type tol = 1e-5;
    constexpr real_type dr  = 1e-5;
    for(int i = 0; i < 1000; ++i)
    {
        for(std::size_t idx=0; idx<3; ++idx)
        {
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dx = uni(mt) * dr;

                mjolnir::math::X(sys.position(idx)) += dx;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::X(sys.position(idx)) += dx;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dx * mjolnir::math::X(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dy = uni(mt) * dr;

                mjolnir::math::Y(sys.position(idx)) += dy;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Y(sys.position(idx)) += dy;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dy * mjolnir::math::Y(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dz = uni(mt) * dr;

                mjolnir::math::Z(sys.position(idx)) += dz;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Z(sys.position(idx)) += dz;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dz * mjolnir::math::Z(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
    }
}
