#define BOOST_TEST_MODULE "test_periodic_gaussian_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/local/PeriodicGaussianPotential.hpp>
#include <mjolnir/math/constants.hpp>

BOOST_AUTO_TEST_CASE(PeriodicGaussian_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi();
    const real_type e  = 2.0;
    const real_type w  = 1.0;
    const real_type r0 = 3.0;

    mjolnir::PeriodicGaussianPotential<real_type> potential(e, w, r0);

    const real_type x_min = -2 * pi;
    const real_type x_max =  2 * pi;
    mjolnir::test::check_potential(potential, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(PeriodicGaussian_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3;
    constexpr real_type tol = 1e-3;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi();
    const real_type e  = 2.0;
    const real_type w  = 1.0;
    const real_type r0 = 3.0;

    mjolnir::PeriodicGaussianPotential<real_type> potential(e, w, r0);

    const real_type x_min = -2 * pi;
    const real_type x_max =  2 * pi;
    mjolnir::test::check_potential(potential, x_min, x_max, tol, h, N);
}
