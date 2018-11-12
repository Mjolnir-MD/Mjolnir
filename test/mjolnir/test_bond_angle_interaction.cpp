#define BOOST_TEST_MODULE "test_bond_angle_interaction"

#include <boost/test/included/unit_test.hpp>
#include <test/util/traits.hpp>
#include <mjolnir/interaction/BondAngleInteraction.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(BondAngleInteraction_force)
{
    typedef mjolnir::test::traits<double> traits;
    constexpr static traits::real_type tol = 1e-7;

    typedef traits::real_type real_type;
    typedef traits::coordinate_type            coord_type;
    typedef traits::boundary_type              boundary_type;
    typedef mjolnir::System<traits>            system_type;
    typedef mjolnir::HarmonicPotential<real_type> harmonic_type;
    typedef mjolnir::BondAngleInteraction<traits, harmonic_type> bond_angle_type;
    typedef bond_angle_type::connection_kind_type connection_kind_type;

    auto normalize = [](const coord_type& v){return v / mjolnir::length(v);};

    const real_type k(1e0);
    const real_type native(mjolnir::math::constants<real_type>::pi * 2.0 / 3.0); // 120 degree

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    const coord_type pos1(1., 0., 0.);
    const coord_type pos2(0., 0., 0.);
    system_type sys(3, boundary_type{});
    sys.at(0) = {1.0, 1.0, pos1,              coord_type(0,0,0), coord_type(0,0,0)};
    sys.at(1) = {1.0, 1.0, pos2,              coord_type(0,0,0), coord_type(0,0,0)};
    sys.at(2) = {1.0, 1.0, coord_type(0,0,0), coord_type(0,0,0), coord_type(0,0,0)};

    const std::size_t N = 1800;
    const real_type dtheta = mjolnir::math::constants<real_type>::pi  / N;
    for(int i = 1; i < N; ++i)
    {
        BOOST_TEST(length(sys[0].position - pos1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(length(sys[1].position - pos2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(length(sys[0].velocity)        == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(length(sys[1].velocity)        == 0.0, boost::test_tools::tolerance(tol));
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
        BOOST_TEST(length(sys[0].position - sys[1].position) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(length(sys[2].position - sys[1].position) == 1.0, boost::test_tools::tolerance(tol));

        const real_type force_strength1 = length(sys[0].force);
        const real_type force_strength3 = length(sys[2].force);
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
        BOOST_TEST(length(sum) == 0.0, boost::test_tools::tolerance(tol));

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_TEST(std::abs(force_strength1) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(std::abs(force_strength3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i < 1200) // narrow
        {
            // perpendicular to radius vector
            const real_type dot1 = dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., -1., 0.);
            BOOST_TEST(length(normalize(sys[0].force) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta + mjolnir::math::constants<real_type>::pi / 2.),
                    sin(theta + mjolnir::math::constants<real_type>::pi / 2.),
                    0.);
            BOOST_TEST(length(normalize(sys[2].force) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i > 1200) // extensive
        {
            const real_type dot1 = dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., 1., 0.);
            BOOST_TEST(length(normalize(sys[0].force) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta - mjolnir::math::constants<real_type>::pi / 2.),
                    sin(theta - mjolnir::math::constants<real_type>::pi / 2.),
                    0.);
            BOOST_TEST(length(normalize(sys[2].force) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }

        // perpendicular to z axis
        BOOST_TEST(sys[0].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[1].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[2].force[2] == 0.0, boost::test_tools::tolerance(tol));
    }
}
