#define BOOST_TEST_MODULE "test_dihedral_angle_interaction"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/potential/HarmonicPotential.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/util/make_unique.hpp>

typedef mjolnir::DefaultTraits traits;
typedef traits::real_type real_type;
typedef traits::coordinate_type coordinate_type;
typedef mjolnir::Particle<coordinate_type> particle_type;
constexpr static traits::real_type tolerance = 1e-5;

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
    typedef traits::real_type real_type;
    typedef traits::coordinate_type coordinate_type;
    typedef mjolnir::Particle<coordinate_type> particle_type;
    typedef mjolnir::ParticleContainer<traits> particle_container_type;
    typedef mjolnir::HarmonicPotential<traits> harmonic_type;
    typedef mjolnir::DihedralAngleInteraction<traits, harmonic_type> dihedral_angle_type;

    const real_type k(1e0);
    const real_type native(M_PI * 2.0 / 3.0);

    typename dihedral_angle_type::container_type potentials;
    std::array<std::size_t, 4> indices{{0,1,2,3}};
    harmonic_type potential(k, native);
    potentials.emplace_back(std::move(indices), std::move(potential));
    dihedral_angle_type inter(std::move(potentials));

    particle_container_type pcon(4);

    const coordinate_type pos1(1e0, 0e0, 1e0);
    const coordinate_type pos2(0e0, 0e0, 1e0);
    const coordinate_type pos3(0e0, 0e0, 0e0);

    pcon[0] = mjolnir::make_particle(1., pos1, zero_vec(), zero_vec());
    pcon[1] = mjolnir::make_particle(1., pos2, zero_vec(), zero_vec());
    pcon[2] = mjolnir::make_particle(1., pos3, zero_vec(), zero_vec());
    pcon[3] = mjolnir::make_particle(1., zero_vec(), zero_vec(), zero_vec());

    const real_type dtheta = M_PI / 1800.0;
    for(int i = -1800; i < 1800; ++i)
    {
        BOOST_CHECK_SMALL(mjolnir::length(pcon[0].position - pos1), tolerance);
        BOOST_CHECK_SMALL(mjolnir::length(pcon[1].position - pos2), tolerance);
        BOOST_CHECK_SMALL(mjolnir::length(pcon[2].position - pos3), tolerance);

        BOOST_CHECK_SMALL(mjolnir::length(pcon[0].velocity), tolerance);
        BOOST_CHECK_SMALL(mjolnir::length(pcon[1].velocity), tolerance);
        BOOST_CHECK_SMALL(mjolnir::length(pcon[2].velocity), tolerance);

        pcon[0].force = zero_vec();
        pcon[1].force = zero_vec();
        pcon[2].force = zero_vec();
        pcon[3].force = zero_vec();

        const real_type theta = i * dtheta;
        const coordinate_type pos4(std::cos(theta), -std::sin(theta), 0e0);
        pcon[3].position = pos4;

        const real_type deriv = potential.derivative(theta);
        const real_type coef = std::abs(deriv);

        inter.calc_force(pcon);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV.
        BOOST_CHECK_CLOSE_FRACTION(mjolnir::length(pcon[1].position - pcon[0].position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(mjolnir::length(pcon[2].position - pcon[1].position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(mjolnir::length(pcon[3].position - pcon[2].position), 1e0, tolerance);

        const real_type force_strength1 = mjolnir::length(pcon[0].force);
        const real_type force_strength3 = mjolnir::length(pcon[3].force);
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
        const coordinate_type sum = pcon[0].force + pcon[1].force + pcon[2].force + pcon[3].force;
        BOOST_CHECK_SMALL(mjolnir::length(sum), tolerance);

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_CHECK_SMALL(mjolnir::length(pcon[0].force), tolerance);
            BOOST_CHECK_SMALL(mjolnir::length(pcon[1].force), tolerance);
            BOOST_CHECK_SMALL(mjolnir::length(pcon[2].force), tolerance);
            BOOST_CHECK_SMALL(mjolnir::length(pcon[3].force), tolerance);
        }
        else
        {
            // perpendicular to radius vector
            const real_type normal1 = dot_product(pcon[0].force, pcon[0].position - pcon[1].position);
            const real_type normal4 = dot_product(pcon[3].force, pcon[2].position - pcon[3].position);
            BOOST_CHECK_SMALL(normal1, tolerance);
            BOOST_CHECK_SMALL(normal4, tolerance);
        }

        // perpendicular to z axis
        BOOST_CHECK_SMALL(pcon[0].force[2], tolerance);
        BOOST_CHECK_SMALL(pcon[1].force[2], tolerance);
        BOOST_CHECK_SMALL(pcon[2].force[2], tolerance);
        BOOST_CHECK_SMALL(pcon[3].force[2], tolerance);
    }
}
