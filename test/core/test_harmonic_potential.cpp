#define BOOST_TEST_MODULE "test_harmonic_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/local/HarmonicPotential.hpp>

BOOST_AUTO_TEST_CASE(Harmonic_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    const real_type k  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::HarmonicPotential<real_type> harmonic(k, r0);

    const real_type x_min = 0.5 * r0;
    const real_type x_max = 1.5 * r0;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = harmonic.potential(x + h);
        const real_type pot2 = harmonic.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = harmonic.derivative(x);

        if(std::abs(pot1 / pot2 - 1.0) < h)
        {
            // pot1 and pot2 are almost the same, thus dpot ~ 0.0.
            // to avoid numerical error in dpot, here it checks `deri ~ 0.0`.
            BOOST_TEST(deri == 0.0, boost::test_tools::tolerance(h));
        }
        else
        {
            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(Harmonic_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3f;

    const real_type k  = 1.0f;
    const real_type r0 = 5.0f;

    mjolnir::HarmonicPotential<real_type> harmonic(k, r0);

    const real_type x_min = 0.5 * r0;
    const real_type x_max = 1.5 * r0;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = harmonic.potential(x + h);
        const real_type pot2 = harmonic.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = harmonic.derivative(x);

        if(std::abs(pot1 / pot2 - 1.0f) < h)
        {
            // pot1 and pot2 are almost the same, thus dpot ~ 0.0.
            // to avoid numerical error in dpot, here it checks `deri ~ 0.0`.
            BOOST_TEST(deri == 0.0f, boost::test_tools::tolerance(h));
        }
        else
        {
            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}
