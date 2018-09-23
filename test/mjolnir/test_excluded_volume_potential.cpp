#define BOOST_TEST_MODULE "test_excluded_volume_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>

BOOST_AUTO_TEST_CASE(EXV_double)
{
    using real_type = double;
    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    mjolnir::ExcludedVolumePotential<real_type, mjolnir::IgnoreNothing>
        exv{epsilon, {sigma, sigma}, {}};

    const real_type x_min = 0.8 * sigma;
    const real_type x_max =
        mjolnir::ExcludedVolumePotential<real_type, mjolnir::IgnoreNothing>::cutoff_ratio * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = exv.potential(0, 1, x + h);
        const real_type pot2 = exv.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = exv.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(EXV_float)
{
    using real_type = float;
    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 0.002;
    constexpr static real_type tol = 0.005;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;

    mjolnir::ExcludedVolumePotential<real_type, mjolnir::IgnoreNothing>
        exv{epsilon, {sigma, sigma}, {}};
    constexpr real_type cutoff =
        mjolnir::ExcludedVolumePotential<real_type, mjolnir::IgnoreNothing>::cutoff_ratio;

    const real_type x_min = 0.8f   * sigma;
    const real_type x_max = cutoff * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = exv.potential(0, 1, x + h);
        const real_type pot2 = exv.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = exv.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

