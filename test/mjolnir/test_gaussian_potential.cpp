#define BOOST_TEST_MODULE "test_gaussian_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>

BOOST_AUTO_TEST_CASE(Gaussian_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    const real_type e  = 2.0;
    const real_type w  = 0.15;
    const real_type r0 = 7.0;

    mjolnir::GaussianPotential<real_type> gaussian(e, w, r0);

    const real_type x_min = 0.5 * r0;
    const real_type x_max = 1.5 * r0;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type pot1 = gaussian.potential(x + h);
        const real_type pot2 = gaussian.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = gaussian.derivative(x);

        if(std::abs(dpot) > h && std::abs(deri) > h)
        {
            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(Gaussian_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3;
    const real_type e  = 2.0;
    const real_type w  = 0.15;
    const real_type r0 = 7.0;

    mjolnir::GaussianPotential<real_type> gaussian(e, w, r0);

    const real_type x_min = 0.5 * r0;
    const real_type x_max = 1.5 * r0;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type pot1 = gaussian.potential(x + h);
        const real_type pot2 = gaussian.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = gaussian.derivative(x);

        if(std::abs(dpot) > h && std::abs(deri) > h)
        {
            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}
