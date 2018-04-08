#define BOOST_TEST_MODULE "test_implicit_membrane_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/external/ImplicitMembranePotential.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(ImplicitMembranePotential_double)
{
    typedef mjolnir::SimulatorTraitsBase<double, mjolnir::UnlimitedBoundary> traits;
    constexpr static std::size_t       N = 10000;
    constexpr static traits::real_type h = 1e-6;
    constexpr static traits::real_type tolerance = 1e-5;

    const traits::real_type thickness = 10.0;
    const traits::real_type interaction_magnitude = 1.0;
    const traits::real_type bend = 1.5;
    const std::vector<traits::real_type> hydrophobicities{1., 0.};

    mjolnir::ImplicitMembranePotential<traits>
        im(thickness, interaction_magnitude, bend, hydrophobicities);

    const traits::real_type cutoff_length = im.max_cutoff_length();
    const traits::real_type z_min = -1 * cutoff_length;
    const traits::real_type z_max = cutoff_length;
    const traits::real_type dz = (z_max - z_min) / N;

    traits::real_type z = z_min;
    for(std::size_t i = 0; i < N; ++i)
    {
        const traits::real_type pot1 = im.potential(0, z + h);
        const traits::real_type pot2 = im.potential(0, z - h);
        const traits::real_type dpot = (pot1 - pot2) / (2 * h);
        const traits::real_type deri = im.derivative(0, z);

        if(std::abs(z) > h)
        {
            if(std::abs(deri) > tolerance)
            {
                BOOST_CHECK_CLOSE_FRACTION(dpot, deri, tolerance);
            }
            else
            {
                BOOST_CHECK_SMALL(deri, tolerance);
            }
        }

        const traits::real_type pot0 = im.potential(1, z);
        const traits::real_type deri0 = im.derivative(1, z);

        BOOST_CHECK_SMALL(pot0, h);
        BOOST_CHECK_SMALL(deri0, h);

        z += dz;
    }
}
