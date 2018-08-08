#define BOOST_TEST_MODULE "test_dihedral_angle_interaction"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/mjolnir/traits.hpp>
#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(DihedralAngle_force)
{
    typedef mjolnir::test::traits<double> traits;
    constexpr static traits::real_type tolerance = 1e-7;

    typedef traits::real_type real_type;
    typedef traits::coordinate_type            coord_type;
    typedef traits::boundary_type              boundary_type;
    typedef mjolnir::System<traits>            system_type;
    typedef system_type::particle_type         particle_type;
    typedef mjolnir::HarmonicPotential<real_type> harmonic_type;
    typedef mjolnir::DihedralAngleInteraction<traits, harmonic_type> dihedral_angle_type;
    typedef dihedral_angle_type::connection_kind_type connection_kind_type;

    auto normalize = [](const coord_type& v){return v / mjolnir::length(v);};

    const real_type k(1e0);
    const real_type native(mjolnir::math::constants<real_type>::pi * 2.0 / 3.0);

    harmonic_type potential{k, native};
    dihedral_angle_type interaction("none", {{ {{0,1,2,3}}, potential}});

    const coord_type pos1(1e0, 0e0, 1e0);
    const coord_type pos2(0e0, 0e0, 1e0);
    const coord_type pos3(0e0, 0e0, 0e0);

    system_type sys(4, boundary_type{});
    sys.at(0) = {1.0, 1.0, pos1,              coord_type(0,0,0), coord_type(0,0,0)};
    sys.at(1) = {1.0, 1.0, pos2,              coord_type(0,0,0), coord_type(0,0,0)};
    sys.at(2) = {1.0, 1.0, pos3,              coord_type(0,0,0), coord_type(0,0,0)};
    sys.at(3) = {1.0, 1.0, coord_type(0,0,0), coord_type(0,0,0), coord_type(0,0,0)};

    const real_type dtheta = mjolnir::math::constants<real_type>::pi / 1800.0;
    for(int i = -1800; i < 1800; ++i)
    {
        BOOST_CHECK_SMALL(mjolnir::length(sys[0].position - pos1), tolerance);
        BOOST_CHECK_SMALL(mjolnir::length(sys[1].position - pos2), tolerance);
        BOOST_CHECK_SMALL(mjolnir::length(sys[2].position - pos3), tolerance);

        BOOST_CHECK_SMALL(mjolnir::length(sys[0].velocity), tolerance);
        BOOST_CHECK_SMALL(mjolnir::length(sys[1].velocity), tolerance);
        BOOST_CHECK_SMALL(mjolnir::length(sys[2].velocity), tolerance);

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
        BOOST_CHECK_CLOSE_FRACTION(mjolnir::length(sys[1].position - sys[0].position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(mjolnir::length(sys[2].position - sys[1].position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(mjolnir::length(sys[3].position - sys[2].position), 1e0, tolerance);

        const real_type force_strength1 = mjolnir::length(sys[0].force);
        const real_type force_strength3 = mjolnir::length(sys[3].force);
        if(i == 1200)
        {
            BOOST_CHECK_SMALL(coef, tolerance);
            BOOST_CHECK_SMALL(force_strength1, tolerance);
            BOOST_CHECK_SMALL(force_strength3, tolerance);
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength1, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength3, tolerance);
        }

        // force applied to center particle is equal to sum of others
        const coord_type sum = sys[0].force + sys[1].force + sys[2].force + sys[3].force;
        BOOST_CHECK_SMALL(mjolnir::length(sum), tolerance);

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_CHECK_SMALL(mjolnir::length(sys[0].force), tolerance);
            BOOST_CHECK_SMALL(mjolnir::length(sys[1].force), tolerance);
            BOOST_CHECK_SMALL(mjolnir::length(sys[2].force), tolerance);
            BOOST_CHECK_SMALL(mjolnir::length(sys[3].force), tolerance);
        }
        else
        {
            // perpendicular to radius vector
            const real_type normal1 = dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type normal4 = dot_product(sys[3].force, sys[2].position - sys[3].position);
            BOOST_CHECK_SMALL(normal1, tolerance);
            BOOST_CHECK_SMALL(normal4, tolerance);
        }

        // perpendicular to z axis
        BOOST_CHECK_SMALL(sys[0].force[2], tolerance);
        BOOST_CHECK_SMALL(sys[1].force[2], tolerance);
        BOOST_CHECK_SMALL(sys[2].force[2], tolerance);
        BOOST_CHECK_SMALL(sys[3].force[2], tolerance);
    }
}
