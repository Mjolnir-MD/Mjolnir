#define BOOST_TEST_MODULE "test_flexible_local_angle_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/FLP/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/math/constants.hpp>

BOOST_AUTO_TEST_CASE(FlexibleLocalAngle_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type k  = 1.0;
    const std::array<real_type, 10> term1{{
        5.00,   1.34,  0.84,   1.17,  0.82,  1.00,  1.27, 1.52,   3.20, 10.00
    }};
    const std::array<real_type, 10> term2{{
        0.00, 151.96, 14.61, -46.89, 39.04, -4.86, -1.86, 8.38, 250.03,  0.00
    }};

    mjolnir::FlexibleLocalAnglePotential<real_type> pot(k, term1, term2);

    const real_type x_min = 1.309;
    const real_type x_max = 2.87979;

    mjolnir::test::check_potential(pot, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(FlexibleLocalAngle_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3;
    constexpr real_type tol = 1e-3;

    const real_type k  = 1.0;
    const std::array<real_type, 10> term1{{
        5.00,   1.34,  0.84,   1.17,  0.82,  1.00,  1.27, 1.52,   3.20, 10.00
    }};
    const std::array<real_type, 10> term2{{
        0.00, 151.96, 14.61, -46.89, 39.04, -4.86, -1.86, 8.38, 250.03,  0.00
    }};

    mjolnir::FlexibleLocalAnglePotential<real_type> pot(k, term1, term2);

    const real_type x_min = 1.309;
    const real_type x_max = 2.87979;

    mjolnir::test::check_potential(pot, x_min, x_max, tol, h, N);
}
