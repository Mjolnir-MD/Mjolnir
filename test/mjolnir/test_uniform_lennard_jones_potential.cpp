#define BOOST_TEST_MODULE "test_uniform_lennard_jones_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(LennardJones_double)
{
    using real_type = double;
    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;

    constexpr real_type sigma   = 3.0;
    constexpr real_type epsilon = 1.0;
    mjolnir::UniformLennardJonesPotential<real_type, mjolnir::IgnoreNothing> lj{
        sigma, epsilon, {}
    };
    constexpr real_type cutoff = mjolnir::UniformLennardJonesPotential<
        real_type, mjolnir::IgnoreNothing>::cutoff_ratio;

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = cutoff * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = lj.potential(0, 1, x + h);
        const real_type pot2 = lj.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = lj.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(LennardJones_float)
{
    using real_type = float;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002;
    constexpr real_type tol = 0.005;

    constexpr real_type sigma   = 3.0;
    constexpr real_type epsilon = 1.0;
    mjolnir::UniformLennardJonesPotential<real_type, mjolnir::IgnoreNothing> lj{
        sigma, epsilon, {}
    };
    constexpr real_type cutoff = mjolnir::UniformLennardJonesPotential<
        real_type, mjolnir::IgnoreNothing>::cutoff_ratio;

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = cutoff * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = lj.potential(0, 1, x + h);
        const real_type pot2 = lj.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = lj.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
    }
}
