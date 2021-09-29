#define BOOST_TEST_MODULE "test_gocontact_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/local/GoContactPotential.hpp>

BOOST_AUTO_TEST_CASE(GoContact_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactPotential<real_type> potential(e, r0);

    const real_type x_min = 0.8 * r0;
    const real_type x_max = 5.0 * r0;

    mjolnir::test::check_potential(potential, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(GoContact_float)
{
    using real_type = double;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3;
    constexpr real_type tol = 1e-3;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactPotential<real_type> potential(e, r0);

    const real_type x_min = 0.8 * r0;
    const real_type x_max = 5.0 * r0;

    mjolnir::test::check_potential(potential, x_min, x_max, tol, h, N);
}
