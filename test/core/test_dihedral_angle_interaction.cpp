#define BOOST_TEST_MODULE "test_dihedral_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/DihedralAngleInteraction.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(DihedralAngle_force)
{
    using traits_type         = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type           = traits_type::real_type;
    using coord_type          = traits_type::coordinate_type;
    using boundary_type       = traits_type::boundary_type;
    using system_type         = mjolnir::System<traits_type>;
    using harmonic_type       = mjolnir::HarmonicPotential<real_type>;
    using dihedral_angle_type = mjolnir::DihedralAngleInteraction<traits_type, harmonic_type>;

    constexpr real_type tol = 1e-7;

    const real_type k(1e0);
    const real_type native(mjolnir::math::constants<real_type>::pi * 2.0 / 3.0);

    harmonic_type potential{k, native};
    dihedral_angle_type interaction("none", {{ {{0,1,2,3}}, potential}});

    const coord_type pos1(1e0, 0e0, 1e0);
    const coord_type pos2(0e0, 0e0, 1e0);
    const coord_type pos3(0e0, 0e0, 0e0);

    system_type sys(4, boundary_type{});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(2).mass = 1.0;
    sys.at(3).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;
    sys.at(2).rmass = 1.0;
    sys.at(3).rmass = 1.0;

    sys.at(0).position = pos1;
    sys.at(1).position = pos2;
    sys.at(2).position = pos3;
    sys.at(3).position = coord_type(0,0,0);
    sys.at(0).velocity = coord_type(0,0,0);
    sys.at(1).velocity = coord_type(0,0,0);
    sys.at(2).velocity = coord_type(0,0,0);
    sys.at(3).velocity = coord_type(0,0,0);
    sys.at(0).force    = coord_type(0,0,0);
    sys.at(1).force    = coord_type(0,0,0);
    sys.at(2).force    = coord_type(0,0,0);
    sys.at(3).force    = coord_type(0,0,0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(2).name  = "X";
    sys.at(3).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";
    sys.at(2).group = "NONE";
    sys.at(3).group = "NONE";

    const real_type dtheta = mjolnir::math::constants<real_type>::pi / 1800.0;
    for(int i = -1800; i < 1800; ++i)
    {
        BOOST_TEST(mjolnir::math::length(sys[0].position - pos1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[1].position - pos2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[2].position - pos3) == 0.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(mjolnir::math::length(sys[0].velocity) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[1].velocity) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[2].velocity) == 0.0, boost::test_tools::tolerance(tol));

        sys[0].force = coord_type(0,0,0);
        sys[1].force = coord_type(0,0,0);
        sys[2].force = coord_type(0,0,0);
        sys[3].force = coord_type(0,0,0);

        const real_type theta = i * dtheta;
        const coord_type pos4(std::cos(theta), -std::sin(theta), 0e0);
        sys[3].position = pos4;

        const real_type deriv = potential.derivative(theta);
        const real_type coef = std::abs(deriv);

        interaction.calc_force(sys);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV.
        BOOST_TEST(mjolnir::math::length(sys[1].position - sys[0].position) == 1e0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[2].position - sys[1].position) == 1e0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[3].position - sys[2].position) == 1e0, boost::test_tools::tolerance(tol));

        const real_type force_strength1 = mjolnir::math::length(sys[0].force);
        const real_type force_strength3 = mjolnir::math::length(sys[3].force);
        if(i == 1200)
        {
            BOOST_TEST(coef            == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(force_strength1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(force_strength3 == 0.0, boost::test_tools::tolerance(tol));
        }
        else
        {
            BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == force_strength3, boost::test_tools::tolerance(tol));
        }

        // force applied to center particle is equal to sum of others
        const coord_type sum = sys[0].force + sys[1].force + sys[2].force + sys[3].force;
        BOOST_TEST(mjolnir::math::length(sum) == 0.0, boost::test_tools::tolerance(tol));

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_TEST(mjolnir::math::length(sys[0].force) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::length(sys[1].force) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::length(sys[2].force) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::length(sys[3].force) == 0.0, boost::test_tools::tolerance(tol));
        }
        else
        {
            // perpendicular to radius vector
            const real_type normal1 = mjolnir::math::dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type normal4 = mjolnir::math::dot_product(sys[3].force, sys[2].position - sys[3].position);
            BOOST_TEST(normal1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(normal4 == 0.0, boost::test_tools::tolerance(tol));
        }

        // perpendicular to z axis
        BOOST_TEST(sys[0].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[1].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[2].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[3].force[2] == 0.0, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(DihedralAngleInteraction_numerical_diff)
{
    using traits_type         = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type           = traits_type::real_type;
    using coord_type          = traits_type::coordinate_type;
    using boundary_type       = traits_type::boundary_type;
    using system_type         = mjolnir::System<traits_type>;
    using harmonic_type       = mjolnir::HarmonicPotential<real_type>;
    using dihedral_angle_type = mjolnir::DihedralAngleInteraction<traits_type, harmonic_type>;

    const real_type k(1e0);
    const real_type native(mjolnir::math::constants<real_type>::pi / 2.0);

    harmonic_type potential{k, native};
    dihedral_angle_type interaction("none", {{ {{0,1,2,3}}, potential}});

    system_type sys(4, boundary_type{});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(2).mass = 1.0;
    sys.at(3).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;
    sys.at(2).rmass = 1.0;
    sys.at(3).rmass = 1.0;

    sys.at(0).position = coord_type(1.0,  0.0, 1.0);
    sys.at(1).position = coord_type(0.0,  0.0, 1.0);
    sys.at(2).position = coord_type(0.0,  0.0, 0.0);
    sys.at(3).position = coord_type(0.0, -1.0, 0.0);
    sys.at(0).velocity = coord_type(0.0,  0.0, 0.0);
    sys.at(1).velocity = coord_type(0.0,  0.0, 0.0);
    sys.at(2).velocity = coord_type(0.0,  0.0, 0.0);
    sys.at(3).velocity = coord_type(0.0,  0.0, 0.0);
    sys.at(0).force    = coord_type(0.0,  0.0, 0.0);
    sys.at(1).force    = coord_type(0.0,  0.0, 0.0);
    sys.at(2).force    = coord_type(0.0,  0.0, 0.0);
    sys.at(3).force    = coord_type(0.0,  0.0, 0.0);

    const auto init = sys;

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(2).name  = "X";
    sys.at(3).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";
    sys.at(2).group = "NONE";
    sys.at(3).group = "NONE";

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    constexpr real_type tol = 1e-5;
    constexpr real_type dr  = 1e-5;
    for(int i = 0; i < 1000; ++i)
    {
        for(std::size_t idx=0; idx<4; ++idx)
        {
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dx = uni(mt) * dr;

                mjolnir::math::X(sys.position(0)) += dx;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::X(sys.position(0)) += dx;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dx * mjolnir::math::X(sys.force(0)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dy = uni(mt) * dr;

                mjolnir::math::Y(sys.position(0)) += dy;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Y(sys.position(0)) += dy;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dy * mjolnir::math::Y(sys.force(0)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dz = uni(mt) * dr;

                mjolnir::math::Z(sys.position(0)) += dz;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Z(sys.position(0)) += dz;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dz * mjolnir::math::Z(sys.force(0)),
                           boost::test_tools::tolerance(tol));
            }
        }
    }
}
