#define BOOST_TEST_MODULE "test_bond_angle_interaction"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/potential/HarmonicPotential.hpp>
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

BOOST_AUTO_TEST_CASE(BondAngleInteraction_force)
{
    const real_type k(1e0);
    const real_type native(M_PI * 2.0 / 3.0); // 120 degree
    std::unique_ptr<mjolnir::LocalPotentialBase<traits>> potential =
        mjolnir::make_unique<mjolnir::HarmonicPotential<traits>>(k, native);
    mjolnir::BondAngleInteraction<traits> inter;

    const coordinate_type pos1(1., 0., 0.);
    const coordinate_type pos2(0., 0., 0.);
    particle_type p1 = mjolnir::make_particle(1., pos1, zero_vec(), zero_vec());
    particle_type p2 = mjolnir::make_particle(1., pos2, zero_vec(), zero_vec());
    particle_type p3 = mjolnir::make_particle(1., zero_vec(), zero_vec(), zero_vec());

    const std::size_t N = 1800;
    const real_type dtheta = M_PI / N;
    for(int i = 1; i < N; ++i)
    {
        BOOST_CHECK_SMALL(length(p1.position - pos1), tolerance);
        BOOST_CHECK_SMALL(length(p2.position - pos2), tolerance);
        BOOST_CHECK_SMALL(length(p1.velocity), tolerance);
        BOOST_CHECK_SMALL(length(p2.velocity), tolerance);
        p1.force = zero_vec();
        p2.force = zero_vec();
        p3.force = zero_vec();

        const real_type theta = i * dtheta;
        const coordinate_type pos3(std::cos(theta), std::sin(theta), 0e0);
        p3.position = pos3;

        const real_type deriv = potential->derivative(theta);
        const real_type coef = std::abs(deriv);

        inter.calc_force(p1, p2, p3, *potential);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV/dtheta.
        BOOST_CHECK_CLOSE_FRACTION(length(p1.position - p2.position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(length(p3.position - p2.position), 1e0, tolerance);

        const real_type force_strength1 = length(p1.force);
        const real_type force_strength3 = length(p3.force);
        if(i == 1200) // most stable point
        {
            BOOST_CHECK_SMALL(coef, tolerance);
            BOOST_CHECK_SMALL(coef, tolerance);
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength1, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength3, tolerance);
        }

        // force applied to center particle is equal to sum of others
        const coordinate_type sum = p1.force + p2.force + p3.force;
        BOOST_CHECK_SMALL(length(sum), tolerance);

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_CHECK_SMALL(std::abs(force_strength1), tolerance);
            BOOST_CHECK_SMALL(std::abs(force_strength3), tolerance);
        }
        else if(i < 1200) // narrow
        {
            // perpendicular to radius vector
            const real_type dot1 = dot_product(p1.force, p1.position - p2.position);
            const real_type dot2 = dot_product(p3.force, p3.position - p2.position);
            BOOST_CHECK_SMALL(dot1, tolerance);
            BOOST_CHECK_SMALL(dot2, tolerance);

            const coordinate_type f1(0., -1., 0.);
            BOOST_CHECK_SMALL(length(normalize(p1.force) - f1), tolerance);

            const coordinate_type f3(cos(theta + M_PI / 2.), sin(theta + M_PI / 2.), 0.);
            BOOST_CHECK_SMALL(length(normalize(p3.force) - f3), tolerance);
        }
        else if(i > 1200) // extensive
        {
            const real_type dot1 = dot_product(p1.force, p1.position - p2.position);
            const real_type dot2 = dot_product(p3.force, p3.position - p2.position);
            BOOST_CHECK_SMALL(dot1, tolerance);
            BOOST_CHECK_SMALL(dot2, tolerance);

            const coordinate_type f1(0., 1., 0.);
            BOOST_CHECK_SMALL(length(normalize(p1.force) - f1), tolerance);

            const coordinate_type f3(cos(theta - M_PI / 2.), sin(theta - M_PI / 2.), 0.);
            BOOST_CHECK_SMALL(length(normalize(p3.force) - f3), tolerance);
        }

        // perpendicular to z axis
        BOOST_CHECK_SMALL(p1.force[2], tolerance);
        BOOST_CHECK_SMALL(p2.force[2], tolerance);
        BOOST_CHECK_SMALL(p3.force[2], tolerance);
    }
}

// BOOST_AUTO_TEST_CASE(BondAngleInteraction_energy)
// {
//     const real_type k(1e0);
//     const real_type native(M_PI * 2.0 / 3.0);
//     LocalPotentialBaseSptr potential =
//         std::make_shared<HarmonicPotential>(k, native);
//
//     InteractionBaseSptr<3> angl(new BondAngleInteraction(potential));
//
//     const coordinate_type pos1(1e0, 0e0, 0e0);
//     const coordinate_type pos2(0e0, 0e0, 0e0);
//     const real_type dtheta = M_PI / 1800.0;
//     for(int i = 0; i < 1000; ++i)
//     {
//         const real_type theta = i * dtheta;
//         const coordinate_type pos3(std::cos(theta), std::sin(theta), 0e0);
//
//         const real_type pot = potential->calc_potential(theta);
//
//         const real_type ene = angl->calc_energy(
//                 std::array<coordinate_type, 3>{{pos1, pos2, pos3}});
//
//         if(i == 1200) // native
//             BOOST_CHECK_SMALL(ene, tolerance);
//         else 
//             BOOST_CHECK_CLOSE_FRACTION(pot, ene, tolerance);
//     }
// }
//

