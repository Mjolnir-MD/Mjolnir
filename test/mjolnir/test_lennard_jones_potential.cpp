#define BOOST_TEST_MODULE "test_lennard_jones_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(LennardJones_double)
{
    using real_type = double;
    constexpr static std::size_t N = 10000;
    constexpr static real_type   h = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    mjolnir::LennardJonesPotential<real_type, mjolnir::IgnoreNothing> lj{
        {{sigma, epsilon}, {sigma, epsilon}}, {}
    };

    const real_type x_min = 0.8 * sigma;
    const real_type x_max =
        mjolnir::LennardJonesPotential<real_type, mjolnir::IgnoreNothing>::cutoff_ratio * sigma;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type pot1 = lj.potential(0, 1, x + h);
        const real_type pot2 = lj.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = lj.derivative(0, 1, x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);

        x += dx;
    }
}

BOOST_AUTO_TEST_CASE(LennardJones_float)
{
    using real_type = float;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002f;
    constexpr real_type tol = 0.005f;

    const real_type sigma   = 3.0f;
    const real_type epsilon = 1.0f;
    mjolnir::LennardJonesPotential<real_type, mjolnir::IgnoreNothing> lj{
        {{sigma, epsilon}, {sigma, epsilon}}, {}
    };
    constexpr real_type cutoff = mjolnir::LennardJonesPotential<
        real_type, mjolnir::IgnoreNothing>::cutoff_ratio;

    const real_type x_min = 0.8f   * sigma;
    const real_type x_max = cutoff * sigma;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type pot1 = lj.potential(0, 1, x + h);
        const real_type pot2 = lj.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = lj.derivative(0, 1, x);

        if(std::abs(deri) > tol)
        {
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, tol);
        }
        else
        {
            BOOST_CHECK_SMALL(deri, tol);
        }
        x += dx;
    }
}


