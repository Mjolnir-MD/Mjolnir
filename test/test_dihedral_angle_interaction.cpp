#define BOOST_TEST_MODULE "test_dihedral_angle_interaction"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/core/HarmonicPotential.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/util/make_unique.hpp>

typedef mjolnir::DefaultTraits traits;
typedef traits::real_type real_type;
typedef traits::coordinate_type coordinate_type;
typedef mjolnir::Particle<coordinate_type> particle_type;
constexpr static traits::real_type tolerance = 1e-7;

traits::coordinate_type zero_vec()
{
    return traits::coordinate_type(0., 0., 0.);
}

traits::coordinate_type normalize(const traits::coordinate_type& vec)
{
    const traits::real_type invl = 1. / mjolnir::length(vec);
    return vec * invl;
}

BOOST_AUTO_TEST_CASE(DihedralAngle_force)
{
    const real_type k(1e0);
    const real_type native(M_PI * 2.0 / 3.0);
    std::unique_ptr<mjolnir::LocalPotentialBase<traits>> potential =
        mjolnir::make_unique<mjolnir::HarmonicPotential<traits>>(k, native);
    mjolnir::DihedralAngleInteraction<traits> inter;

    const coordinate_type pos1(1e0, 0e0, 1e0);
    const coordinate_type pos2(0e0, 0e0, 1e0);
    const coordinate_type pos3(0e0, 0e0, 0e0);

    particle_type p1 = mjolnir::make_particle(1., pos1, zero_vec(), zero_vec());
    particle_type p2 = mjolnir::make_particle(1., pos2, zero_vec(), zero_vec());
    particle_type p3 = mjolnir::make_particle(1., pos3, zero_vec(), zero_vec());

    const real_type dtheta = M_PI / 1800.0;
    for(int i = -1800; i < 1800; ++i)
    {
        BOOST_CHECK_SMALL(length(p1.position - pos1), tolerance);
        BOOST_CHECK_SMALL(length(p2.position - pos2), tolerance);
        BOOST_CHECK_SMALL(length(p3.position - pos3), tolerance);

        BOOST_CHECK_SMALL(length(p1.velocity), tolerance);
        BOOST_CHECK_SMALL(length(p2.velocity), tolerance);
        BOOST_CHECK_SMALL(length(p3.velocity), tolerance);

        p1.force = zero_vec();
        p2.force = zero_vec();
        p3.force = zero_vec();

        const real_type theta = i * dtheta;
        const coordinate_type pos4(std::cos(theta), -std::sin(theta), 0e0);
        particle_type p4 = mjolnir::make_particle(1., pos4, zero_vec(), zero_vec());

        const real_type deriv = potential->derivative(theta);
        const real_type coef = std::abs(deriv);

        inter.calc_force(p1, p2, p3, p4, *potential);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV.
        BOOST_CHECK_CLOSE_FRACTION(length(p2.position - p1.position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(length(p3.position - p2.position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(length(p4.position - p3.position), 1e0, tolerance);

        const real_type force_strength1 = length(p1.force);
        const real_type force_strength3 = length(p4.force);
        BOOST_CHECK_CLOSE_FRACTION(coef, force_strength1, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(coef, force_strength3, tolerance);

        // force applied to center particle is equal to sum of others
        const coordinate_type sum = p1.force + p2.force + p3.force + p4.force;
        BOOST_CHECK_SMALL(length(sum), tolerance);

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_CHECK_SMALL(length(p1.force), tolerance);
            BOOST_CHECK_SMALL(length(p2.force), tolerance);
            BOOST_CHECK_SMALL(length(p3.force), tolerance);
            BOOST_CHECK_SMALL(length(p4.force), tolerance);
        }
        else
        {
            // perpendicular to radius vector
            const real_type normal1 = dot_product(p1.force, p1.position - p2.position);
            const real_type normal4 = dot_product(p4.force, p3.position - p4.position);
            BOOST_CHECK_SMALL(normal1, tolerance);
            BOOST_CHECK_SMALL(normal4, tolerance);
        }

        // perpendicular to z axis
        BOOST_CHECK_SMALL(p1.force[2], tolerance);
        BOOST_CHECK_SMALL(p2.force[2], tolerance);
        BOOST_CHECK_SMALL(p3.force[2], tolerance);
        BOOST_CHECK_SMALL(p4.force[2], tolerance);
    }
}
