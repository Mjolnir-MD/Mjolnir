#define BOOST_TEST_MODULE "test_fast_inv_sqrt"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/math/math.hpp>

#include <random>
#include <chrono>
#include <array>
#include <iostream>

constexpr static std::size_t N = 10000;

BOOST_AUTO_TEST_CASE(rsqrt_double)
{
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(0.0001, 100.0);
    const double tolerance = 1.0 * std::pow(2.0, -14);
    for(std::size_t i=0; i < N; ++i)
    {
        const double x = uni(mt);
        const double x2 = x * x;
        const double x_inv = 1.0 / x;
        BOOST_CHECK_CLOSE_FRACTION(x_inv, mjolnir::rsqrt(x2), tolerance);
    }
}

BOOST_AUTO_TEST_CASE(rsqrt_float)
{
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<float> uni(0.0001, 100.0);
    const float tolerance = 1.5 * std::pow(2.0, -12);
    for(std::size_t i=0; i < N; ++i)
    {
        const float x = uni(mt);
        const float x2 = x * x;
        const float x_inv = 1.0f / x;
        BOOST_CHECK_CLOSE_FRACTION(x_inv, mjolnir::rsqrt(x2), tolerance);
    }
}
