#define BOOST_TEST_MODULE "test_potential_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/local/CosinePotential.hpp>
#include <mjolnir/math/constants.hpp>

BOOST_AUTO_TEST_CASE(CosinePotential_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi;
    const real_type    k = 2.0;
    const std::int32_t n = 2;
    const real_type   v0 = 3.0;

    mjolnir::CosinePotential<real_type> potential(k, n, v0);

    const real_type x_min = -pi;
    const real_type x_max =  pi;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type pot1 = potential.potential(x + h);
        const real_type pot2 = potential.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = potential.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(CosinePotential_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi;
    const real_type    k = 2.0;
    const std::int32_t n = 2;
    const real_type   v0 = 3.0;

    mjolnir::CosinePotential<real_type> potential(k, n, v0);

    const real_type x_min = -2 * pi;
    const real_type x_max =  2 * pi;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type pot1 = potential.potential(x + h);
        const real_type pot2 = potential.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = potential.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}
