#define BOOST_TEST_MODULE "test_dihedral_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/DihedralAngleInteraction.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/forcefield/local/ClementiDihedralPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(DihedralAngleInteraction_numerical_diff)
{
    namespace test = mjolnir::test;
    mjolnir::LoggerManager::set_default_logger("test_dihedral_interaction.log");

    using traits_type         = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type           = traits_type::real_type;
    using coord_type          = traits_type::coordinate_type;
    using boundary_type       = traits_type::boundary_type;
    using system_type         = mjolnir::System<traits_type>;
    using potential_type       = mjolnir::ClementiDihedralPotential<real_type>;
    using dihedral_angle_type = mjolnir::DihedralAngleInteraction<traits_type, potential_type>;

    constexpr auto pi = mjolnir::math::constants<real_type>::pi();

    const real_type k1(1e0);
    const real_type k3(1e0);
    const real_type native(pi / 2.0);

    std::mt19937 mt(123456789);

    potential_type potential{k1, k3, native};
    dihedral_angle_type interaction("none", {{ {{0,1,2,3}}, potential}});

    for(std::size_t i=0; i<1000; ++i)
    {
        const real_type theta = i * 0.001 * (2.0 * pi);

        system_type sys(4, boundary_type{});
        test::clear_everything(sys);

        sys.position(0) = coord_type( 2.0,  std::cos(theta), std::sin(theta));
        sys.position(1) = coord_type( 1.0,  0.0, 0.0);
        sys.position(2) = coord_type( 0.0,  0.0, 0.0);
        sys.position(3) = coord_type(-1.0, -1.0, 0.0);

        test::apply_random_rotation(sys, mt);
        test::apply_random_perturbation(sys, mt, 0.01);

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, interaction, tol, dr);
        test::check_virial(sys, interaction, tol);
        test::check_force_and_energy(sys, interaction, tol);
    }
}

BOOST_AUTO_TEST_CASE(DihedralAngle_force)
{
    mjolnir::LoggerManager::set_default_logger("test_dihedral_interaction.log");

    using traits_type         = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type           = traits_type::real_type;
    using coord_type          = traits_type::coordinate_type;
    using boundary_type       = traits_type::boundary_type;
    using system_type         = mjolnir::System<traits_type>;
    using potential_type      = mjolnir::ClementiDihedralPotential<real_type>;
    using dihedral_angle_type = mjolnir::DihedralAngleInteraction<traits_type, potential_type>;

    constexpr real_type tol = 1e-7;

    const real_type k1(1e0);
    const real_type k3(1e0);
    const real_type native(mjolnir::math::constants<real_type>::pi() * 2.0 / 3.0);

    potential_type potential{k1, k3, native};
    dihedral_angle_type interaction("none", {{ {{0,1,2,3}}, potential}});

    const coord_type pos1(1e0, 0e0, 1e0);
    const coord_type pos2(0e0, 0e0, 1e0);
    const coord_type pos3(0e0, 0e0, 0e0);

    system_type sys(4, boundary_type{});

    sys.mass(0) = 1.0;
    sys.mass(1) = 1.0;
    sys.mass(2) = 1.0;
    sys.mass(3) = 1.0;
    sys.rmass(0) = 1.0;
    sys.rmass(1) = 1.0;
    sys.rmass(2) = 1.0;
    sys.rmass(3) = 1.0;

    sys.position(0) = pos1;
    sys.position(1) = pos2;
    sys.position(2) = pos3;
    sys.position(3) = coord_type(0,0,0);
    sys.velocity(0) = coord_type(0,0,0);
    sys.velocity(1) = coord_type(0,0,0);
    sys.velocity(2) = coord_type(0,0,0);
    sys.velocity(3) = coord_type(0,0,0);
    sys.force(0)    = coord_type(0,0,0);
    sys.force(1)    = coord_type(0,0,0);
    sys.force(2)    = coord_type(0,0,0);
    sys.force(3)    = coord_type(0,0,0);

    sys.name(0)  = "X";
    sys.name(1)  = "X";
    sys.name(2)  = "X";
    sys.name(3)  = "X";
    sys.group(0) = "NONE";
    sys.group(1) = "NONE";
    sys.group(2) = "NONE";
    sys.group(3) = "NONE";

    const real_type dtheta = mjolnir::math::constants<real_type>::pi() / 1800.0;
    for(int i = -1800; i < 1800; ++i)
    {
        BOOST_TEST(mjolnir::math::length(sys.position(0) - pos1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.position(1) - pos2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.position(2) - pos3) == 0.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(mjolnir::math::length(sys.velocity(0)) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.velocity(1)) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.velocity(2)) == 0.0, boost::test_tools::tolerance(tol));

        sys.force(0) = coord_type(0,0,0);
        sys.force(1) = coord_type(0,0,0);
        sys.force(2) = coord_type(0,0,0);
        sys.force(3) = coord_type(0,0,0);

        const real_type theta = i * dtheta;
        const coord_type pos4(std::cos(theta), -std::sin(theta), 0e0);
        sys.position(3) = pos4;

        const real_type deriv = potential.derivative(theta);
        const real_type coef = std::abs(deriv);

        interaction.calc_force(sys);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV.
        BOOST_TEST(mjolnir::math::length(sys.position(1) - sys.position(0)) == 1e0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.position(2) - sys.position(1)) == 1e0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.position(3) - sys.position(2)) == 1e0, boost::test_tools::tolerance(tol));

        const real_type force_strength1 = mjolnir::math::length(sys.force(0));
        const real_type force_strength3 = mjolnir::math::length(sys.force(3));
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
        const coord_type sum = sys.force(0) + sys.force(1) + sys.force(2) + sys.force(3);
        BOOST_TEST(mjolnir::math::length(sum) == 0.0, boost::test_tools::tolerance(tol));

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_TEST(mjolnir::math::length(sys.force(0)) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::length(sys.force(1)) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::length(sys.force(2)) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::length(sys.force(3)) == 0.0, boost::test_tools::tolerance(tol));
        }
        else
        {
            // perpendicular to radius vector
            const real_type normal1 = mjolnir::math::dot_product(sys.force(0), sys.position(0) - sys.position(1));
            const real_type normal4 = mjolnir::math::dot_product(sys.force(3), sys.position(2) - sys.position(3));
            BOOST_TEST(normal1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(normal4 == 0.0, boost::test_tools::tolerance(tol));
        }

        // perpendicular to z axis
        BOOST_TEST(sys.force(0)[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(1)[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(2)[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(3)[2] == 0.0, boost::test_tools::tolerance(tol));
    }
}
