#define BOOST_TEST_MODULE "test_harmonic_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>

BOOST_AUTO_TEST_CASE(Harmonic_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type k  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::HarmonicPotential<real_type> harmonic(k, r0);

    const real_type x_min = 0.5 * r0;
    const real_type x_max = 1.5 * r0;

    mjolnir::test::check_potential(harmonic, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(Harmonic_float)
{
    namespace test = mjolnir::test;
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3f;
    constexpr real_type tol = 1e-3f;

    const real_type k  = 1.0f;
    const real_type r0 = 5.0f;

    mjolnir::HarmonicPotential<real_type> harmonic(k, r0);

    const real_type x_min = 0.5 * r0;
    const real_type x_max = 1.5 * r0;

    mjolnir::test::check_potential(harmonic, x_min, x_max, tol, h, N);
}
