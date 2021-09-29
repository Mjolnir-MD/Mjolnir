#define BOOST_TEST_MODULE "test_bond_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/BondAngleInteraction.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(BondAngleInteraction_numerical_diff)
{
    namespace test = mjolnir::test;

    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = traits_type::real_type;
    using coord_type      = traits_type::coordinate_type;
    using boundary_type   = traits_type::boundary_type;
    using system_type     = mjolnir::System<traits_type>;
    using harmonic_type   = mjolnir::HarmonicPotential<real_type>;
    using bond_angle_type = mjolnir::BondAngleInteraction<traits_type, harmonic_type>;

    const real_type k(10.0);
    const real_type native(mjolnir::math::constants<real_type>::pi() / 3.0);

    std::mt19937 mt(123456789);

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    for(int i = 0; i < 1000; ++i)
    {
        system_type sys(3, boundary_type{});
        test::clear_everything(sys);

        sys.position(0) = coord_type(1.0, 0.0, 0.0);
        sys.position(1) = coord_type(0.0, 0.0, 1.0);
        sys.position(2) = coord_type(1.0, 1.0, 1.0);

        test::apply_random_rotation(sys, mt);
        test::apply_random_perturbation(sys, mt, 0.01);

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, interaction, tol, dr);
        test::check_virial(sys, interaction, tol);
        test::check_force_and_virial(sys, interaction, tol);
        test::check_force_and_energy(sys, interaction, tol);
        test::check_force_energy_virial(sys, interaction, tol);
    }
}

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
    const real_type native(mjolnir::math::constants<real_type>::pi() * 2.0 / 3.0); // 120 degree

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    const coord_type pos1(1., 0., 0.);
    const coord_type pos2(0., 0., 0.);
    system_type sys(3, boundary_type{});

    sys.mass(0) = 1.0;
    sys.mass(1) = 1.0;
    sys.mass(2) = 1.0;
    sys.rmass(0) = 1.0;
    sys.rmass(1) = 1.0;
    sys.rmass(2) = 1.0;

    sys.position(0) = pos1;
    sys.position(1) = pos2;
    sys.position(2) = coord_type(0,0,0);
    sys.velocity(0) = coord_type(0,0,0);
    sys.velocity(1) = coord_type(0,0,0);
    sys.velocity(2) = coord_type(0,0,0);
    sys.force(0)    = coord_type(0,0,0);
    sys.force(1)    = coord_type(0,0,0);
    sys.force(2)    = coord_type(0,0,0);

    sys.name(0)  = "X";
    sys.name(1)  = "X";
    sys.name(2)  = "X";
    sys.group(0) = "NONE";
    sys.group(1) = "NONE";
    sys.group(2) = "NONE";

    const int N = 1800;
    const real_type dtheta = mjolnir::math::constants<real_type>::pi()  / N;
    for(int i = 1; i < N; ++i)
    {
        BOOST_TEST(mjolnir::math::length(sys.position(0) - pos1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.position(1) - pos2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.velocity(0))        == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.velocity(1))        == 0.0, boost::test_tools::tolerance(tol));
        sys.force(0) = coord_type(0,0,0);
        sys.force(1) = coord_type(0,0,0);
        sys.force(2) = coord_type(0,0,0);

        const real_type theta = i * dtheta;
        const coord_type pos3(std::cos(theta), std::sin(theta), 0e0);
        sys.position(2) = pos3;

        const real_type deriv = potential.derivative(theta);
        const real_type coef = std::abs(deriv);

        interaction.calc_force(sys);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV/dtheta.
        BOOST_TEST(mjolnir::math::length(sys.position(0) - sys.position(1)) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys.position(2) - sys.position(1)) == 1.0, boost::test_tools::tolerance(tol));

        const real_type force_strength1 = mjolnir::math::length(sys.force(0));
        const real_type force_strength3 = mjolnir::math::length(sys.force(2));
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

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_TEST(std::abs(force_strength1) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(std::abs(force_strength3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i < 1200) // narrow
        {
            // perpendicular to radius vector
            const real_type dot1 = mjolnir::math::dot_product(sys.force(0), sys.position(0) - sys.position(1));
            const real_type dot2 = mjolnir::math::dot_product(sys.force(2), sys.position(2) - sys.position(1));
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., -1., 0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys.force(0)) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta + mjolnir::math::constants<real_type>::pi() / 2.),
                    sin(theta + mjolnir::math::constants<real_type>::pi() / 2.),
                    0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys.force(2)) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i > 1200) // extensive
        {
            const real_type dot1 = mjolnir::math::dot_product(sys.force(0), sys.position(0) - sys.position(1));
            const real_type dot2 = mjolnir::math::dot_product(sys.force(2), sys.position(2) - sys.position(1));
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., 1., 0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys.force(0)) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta - mjolnir::math::constants<real_type>::pi() / 2.),
                    sin(theta - mjolnir::math::constants<real_type>::pi() / 2.),
                    0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys.force(2)) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }

        // perpendicular to z axis
        BOOST_TEST(mjolnir::math::Z(sys.force(0)) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(sys.force(1)) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(sys.force(2)) == 0.0, boost::test_tools::tolerance(tol));
    }
}
