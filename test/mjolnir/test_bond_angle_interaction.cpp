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
    typedef traits::real_type real_type;
    typedef traits::coordinate_type coordinate_type;
    typedef mjolnir::Particle<coordinate_type> particle_type;
    typedef mjolnir::ParticleContainer<traits> particle_container_type;
    typedef mjolnir::HarmonicPotential<traits> harmonic_type;
    typedef mjolnir::BondAngleInteraction<traits, harmonic_type> bond_angle_type;

    const real_type k(1e0);
    const real_type native(M_PI * 2.0 / 3.0); // 120 degree

    typename bond_angle_type::container_type potentials;
    std::array<std::size_t, 3> indices{{0,1,2}};
    harmonic_type potential(k, native);
    potentials.emplace_back(std::move(indices), std::move(potential));
    bond_angle_type inter(std::move(potentials));

    const coordinate_type pos1(1., 0., 0.);
    const coordinate_type pos2(0., 0., 0.);
    particle_container_type pcon(3);
    pcon[0] = mjolnir::make_particle(1., pos1,       zero_vec(), zero_vec());
    pcon[1] = mjolnir::make_particle(1., pos2,       zero_vec(), zero_vec());
    pcon[2] = mjolnir::make_particle(1., zero_vec(), zero_vec(), zero_vec());

    const std::size_t N = 1800;
    const real_type dtheta = M_PI / N;
    for(int i = 1; i < N; ++i)
    {
        BOOST_CHECK_SMALL(length(pcon[0].position - pos1), tolerance);
        BOOST_CHECK_SMALL(length(pcon[1].position - pos2), tolerance);
        BOOST_CHECK_SMALL(length(pcon[0].velocity), tolerance);
        BOOST_CHECK_SMALL(length(pcon[1].velocity), tolerance);
        pcon[0].force = zero_vec();
        pcon[1].force = zero_vec();
        pcon[2].force = zero_vec();

        const real_type theta = i * dtheta;
        const coordinate_type pos3(std::cos(theta), std::sin(theta), 0e0);
        pcon[2].position = pos3;

        const real_type deriv = potential.derivative(theta);
        const real_type coef = std::abs(deriv);

        inter.calc_force(pcon);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV/dtheta.
        BOOST_CHECK_CLOSE_FRACTION(length(pcon[0].position - pcon[1].position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(length(pcon[2].position - pcon[1].position), 1e0, tolerance);

        const real_type force_strength1 = length(pcon[0].force);
        const real_type force_strength3 = length(pcon[2].force);
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
        const coordinate_type sum = pcon[0].force + pcon[1].force + pcon[2].force;
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
            const real_type dot1 = dot_product(pcon[0].force, pcon[0].position - pcon[1].position);
            const real_type dot2 = dot_product(pcon[2].force, pcon[2].position - pcon[1].position);
            BOOST_CHECK_SMALL(dot1, tolerance);
            BOOST_CHECK_SMALL(dot2, tolerance);

            const coordinate_type f1(0., -1., 0.);
            BOOST_CHECK_SMALL(length(normalize(pcon[0].force) - f1), tolerance);

            const coordinate_type f3(cos(theta + M_PI / 2.), sin(theta + M_PI / 2.), 0.);
            BOOST_CHECK_SMALL(length(normalize(pcon[2].force) - f3), tolerance);
        }
        else if(i > 1200) // extensive
        {
            const real_type dot1 = dot_product(pcon[0].force, pcon[0].position - pcon[1].position);
            const real_type dot2 = dot_product(pcon[2].force, pcon[2].position - pcon[1].position);
            BOOST_CHECK_SMALL(dot1, tolerance);
            BOOST_CHECK_SMALL(dot2, tolerance);

            const coordinate_type f1(0., 1., 0.);
            BOOST_CHECK_SMALL(length(normalize(pcon[0].force) - f1), tolerance);

            const coordinate_type f3(cos(theta - M_PI / 2.), sin(theta - M_PI / 2.), 0.);
            BOOST_CHECK_SMALL(length(normalize(pcon[2].force) - f3), tolerance);
        }

        // perpendicular to z axis
        BOOST_CHECK_SMALL(pcon[0].force[2], tolerance);
        BOOST_CHECK_SMALL(pcon[1].force[2], tolerance);
        BOOST_CHECK_SMALL(pcon[2].force[2], tolerance);
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

