#define BOOST_TEST_MODULE "test_harmonic_restraint_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/external/HarmonicRestraintPotential.hpp>

BOOST_AUTO_TEST_CASE(Harmonic_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    mjolnir::HarmonicRestraintPotential<real_type> harmonic({
            {1.0, 0.0}, // k = 1.0, r0 = 0.0
            {1.0, 1.0}, // k = 1.0, r0 = 1.0
            {2.0, 0.0}, // k = 2.0, r0 = 0.0
            });

    for(std::size_t j=0; j<3; ++j)
    {
        const real_type r0 = harmonic.parameters()[j].second;

        const real_type x_min = (r0 != 0.0) ? 0.5 * r0 : h;
        const real_type x_max = (r0 != 0.0) ? 1.5 * r0 : 1.0;
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x    = x_min + i * dx;
            const real_type pot1 = harmonic.potential(j, x + h);
            const real_type pot2 = harmonic.potential(j, x - h);
            const real_type dpot = (pot1 - pot2) / (2 * h);
            const real_type deri = harmonic.derivative(j, x);

            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(Harmonic_float)
{
    using real_type = float;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-3f;

    mjolnir::HarmonicRestraintPotential<real_type> harmonic({
            {1.0f, 0.0f}, // k = 1.0, r0 = 0.0
            {1.0f, 1.0f}, // k = 1.0, r0 = 1.0
            {2.0f, 0.0f}, // k = 2.0, r0 = 0.0
            });

    for(std::size_t j=0; j<3; ++j)
    {
        const real_type r0 = harmonic.parameters()[j].second;

        const real_type x_min = (r0 != 0.0f) ? 0.5f * r0 : h;
        const real_type x_max = (r0 != 0.0f) ? 1.5f * r0 : 1.0f;
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x    = x_min + i * dx;
            const real_type pot1 = harmonic.potential(j, x + h);
            const real_type pot2 = harmonic.potential(j, x - h);
            const real_type dpot = (pot1 - pot2) / (2 * h);
            const real_type deri = harmonic.derivative(j, x);

            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}
