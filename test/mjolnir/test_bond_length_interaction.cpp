#define BOOST_TEST_MODULE "test_bond_length_interaction"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/potential/HarmonicPotential.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/util/make_unique.hpp>

typedef mjolnir::DefaultTraits traits;
constexpr static traits::real_type tolerance = 1e-8;

traits::coordinate_type zero_vec()
{
    return traits::coordinate_type(0., 0., 0.);
}

traits::coordinate_type normalize(const traits::coordinate_type& vec)
{
    const traits::real_type invl = 1. / mjolnir::length(vec);
    return vec * invl;
}

BOOST_AUTO_TEST_CASE(BondLength_calc_force)
{
    typedef traits::real_type real_type;
    typedef traits::coordinate_type coordinate_type;
    typedef mjolnir::Particle<coordinate_type> particle_type;

    const real_type k(100.);
    const real_type native(2.0);
    std::unique_ptr<mjolnir::LocalPotentialBase<traits>> potential =
        mjolnir::make_unique<mjolnir::HarmonicPotential<traits>>(k, native);
    mjolnir::BondLengthInteraction<traits> inter;

    particle_type p1 = mjolnir::make_particle(1., zero_vec(), zero_vec(), zero_vec());
    particle_type p2 = mjolnir::make_particle(1., zero_vec(), zero_vec(), zero_vec());
    std::array<particle_type*, 2> ps;
    ps[0] = &p1;
    ps[1] = &p2;

    const real_type dr = 1e-3;
    real_type dist = 1e0;
    for(int i = 0; i < 2000; ++i)
    {
        p1.position = zero_vec();
        p2.position = zero_vec();

        p1.force = zero_vec();
        p2.force = zero_vec();

        p2.position[0] = dist;

        const real_type deriv = potential->derivative(dist);// dV/dr: f=-dV/dr
        const real_type coef  = std::abs(deriv);

        inter.calc_force(ps, *potential);

        const real_type force_strength1 = mjolnir::length(p1.force);
        const real_type force_strength2 = mjolnir::length(p2.force);

        BOOST_CHECK_CLOSE_FRACTION(coef, force_strength1, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(coef, force_strength2, tolerance);

        // direction
        if(i == 1000) // most stable point
        {
            BOOST_CHECK_SMALL(force_strength1, tolerance);
            BOOST_CHECK_SMALL(force_strength2, tolerance);
        }
        else if(i < 1000) // repulsive
        {
            const real_type dir1 =
                mjolnir::dot_product(normalize(p1.force),
                                     normalize(p1.position - p2.position));
            const real_type dir2 =
                mjolnir::dot_product(normalize(p2.force),
                                     normalize(p2.position - p1.position));

            BOOST_CHECK_CLOSE_FRACTION(dir1, 1e0, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(dir2, 1e0, tolerance);
        }
        else if(i > 1000) // attractive
        {
            const real_type dir1 =
                mjolnir::dot_product(normalize(p1.force),
                                     normalize(p2.position - p1.position));
            const real_type dir2 =
                mjolnir::dot_product(normalize(p2.force),
                                     normalize(p1.position - p2.position));

            BOOST_CHECK_CLOSE_FRACTION(dir1, 1e0, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(dir2, 1e0, tolerance);
        }

        BOOST_CHECK_SMALL(length(p1.force + p2.force), tolerance);

        dist += dr;
    }
}
