#define BOOST_TEST_MODULE "test_bond_length_interaction"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(BondLength_calc_force)
{
    typedef mjolnir::SimulatorTraitsBase<double, mjolnir::UnlimitedBoundary> traits;
    constexpr static traits::real_type tolerance = 1e-8;

    typedef traits::real_type real_type;
    typedef traits::coordinate_type            coord_type;
    typedef traits::boundary_type              boundary_type;
    typedef traits::particle_type              particle_type;
    typedef mjolnir::System<traits>            system_type;
    typedef mjolnir::HarmonicPotential<traits> harmonic_type;
    typedef mjolnir::BondLengthInteraction<traits, harmonic_type> bond_length_type;
    typedef bond_length_type::connection_kind_type connection_kind_type;

    auto normalize = [](const coord_type& v){return v / mjolnir::length(v);};

    const real_type k(100.);
    const real_type native(2.0);

    harmonic_type    potential(k, native);
    bond_length_type interaction(connection_kind_type::none, {{ {{0,1}}, potential}});

    std::vector<particle_type> ps{
        {1., coord_type(0,0,0), coord_type(0,0,0), coord_type(0,0,0)},
        {1., coord_type(0,0,0), coord_type(0,0,0), coord_type(0,0,0)}
    };
    system_type sys(std::move(ps), boundary_type{});

    const real_type dr = 1e-3;
    real_type dist = 1e0;
    for(int i = 0; i < 2000; ++i)
    {
        sys[0].position = coord_type(0,0,0);
        sys[1].position = coord_type(0,0,0);
        sys[0].force    = coord_type(0,0,0);
        sys[1].force    = coord_type(0,0,0);
        sys[1].position[0] = dist;

        const real_type deriv = potential.derivative(dist);
        const real_type coef  = std::abs(deriv);

        interaction.calc_force(sys);

        const real_type force_strength1 = mjolnir::length(sys[0].force);
        const real_type force_strength2 = mjolnir::length(sys[1].force);


        // direction
        if(i == 1000) // most stable point
        {
            BOOST_CHECK_SMALL(force_strength1, tolerance);
            BOOST_CHECK_SMALL(force_strength2, tolerance);
        }
        else if(i < 1000) // repulsive
        {
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength1, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength2, tolerance);

            const real_type dir1 =
                mjolnir::dot_product(normalize(sys[0].force),
                                     normalize(sys[0].position - sys[1].position));
            const real_type dir2 =
                mjolnir::dot_product(normalize(sys[1].force),
                                     normalize(sys[1].position - sys[0].position));

            BOOST_CHECK_CLOSE_FRACTION(dir1, 1e0, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(dir2, 1e0, tolerance);
        }
        else if(i > 1000) // attractive
        {
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength1, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength2, tolerance);

            const real_type dir1 =
                mjolnir::dot_product(normalize(sys[0].force),
                                     normalize(sys[1].position - sys[0].position));
            const real_type dir2 =
                mjolnir::dot_product(normalize(sys[1].force),
                                     normalize(sys[0].position - sys[1].position));

            BOOST_CHECK_CLOSE_FRACTION(dir1, 1e0, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(dir2, 1e0, tolerance);
        }

        BOOST_CHECK_SMALL(length(sys[0].force + sys[1].force), tolerance);

        dist += dr;
    }
}
