#define BOOST_TEST_MODULE "test_periodic_gaussian_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/potential/local/PeriodicGaussianPotential.hpp>
#include <mjolnir/math/constants.hpp>

BOOST_AUTO_TEST_CASE(PeriodicGaussian_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi;
    const real_type e  = 2.0;
    const real_type w  = 1.0;
    const real_type r0 = 3.0;

    mjolnir::PeriodicGaussianPotential<real_type> gaussian(e, w, r0);

    const real_type x_min = -2 * pi;
    const real_type x_max =  2 * pi;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type pot1 = gaussian.potential(x + h);
        const real_type pot2 = gaussian.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = gaussian.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(PeriodicGaussian_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi;
    const real_type e  = 2.0;
    const real_type w  = 1.0;
    const real_type r0 = 3.0;

    mjolnir::PeriodicGaussianPotential<real_type> gaussian(e, w, r0);

    const real_type x_min = -2 * pi;
    const real_type x_max =  2 * pi;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type pot1 = gaussian.potential(x + h);
        const real_type pot2 = gaussian.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = gaussian.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}
