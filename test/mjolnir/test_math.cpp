#define BOOST_TEST_MODULE "test_math_functions"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <boost/mpl/list.hpp>
#include <mjolnir/math/math.hpp>

#include <random>
#include <chrono>
#include <array>
#include <iostream>

typedef boost::mpl::list<double, float> test_targets;
constexpr std::size_t N = 10000;

static_assert(mjolnir::compiletime::pow( 2.0, 4) == 16.0, "mjolnir::compiletime::pow");
static_assert(mjolnir::compiletime::pow(-2.0, 4) == 16.0, "mjolnir::compiletime::pow");
static_assert(mjolnir::compiletime::pow(-2.0, 3) == -8.0, "mjolnir::compiletime::pow");

static_assert(mjolnir::compiletime::abs( 1.0) == 1.0, "mjolnir::compiletime::abs");
static_assert(mjolnir::compiletime::abs(-1.0) == 1.0, "mjolnir::compiletime::abs");
static_assert(mjolnir::compiletime::min( 1.0, 2.0) ==  1.0, "mjolnir::compiletime::min");
static_assert(mjolnir::compiletime::min(-1.0, 0.1) == -1.0, "mjolnir::compiletime::min");
static_assert(mjolnir::compiletime::min(2.0,  1.0) ==  1.0, "mjolnir::compiletime::min");
static_assert(mjolnir::compiletime::min(0.1, -1.0) == -1.0, "mjolnir::compiletime::min");
static_assert(mjolnir::compiletime::max( 1.0, 2.0) ==  2.0, "mjolnir::compiletime::max");
static_assert(mjolnir::compiletime::max(-1.0, 0.1) ==  0.1, "mjolnir::compiletime::max");
static_assert(mjolnir::compiletime::max(2.0,  1.0) ==  2.0, "mjolnir::compiletime::max");
static_assert(mjolnir::compiletime::max(0.1, -1.0) ==  0.1, "mjolnir::compiletime::max");

static_assert(mjolnir::compiletime::clamp( 0.1, -1.0, 1.0) ==  0.1, "mjolnir::compiletime::clamp");
static_assert(mjolnir::compiletime::clamp( 2.0, -1.0, 1.0) ==  1.0, "mjolnir::compiletime::clamp");
static_assert(mjolnir::compiletime::clamp(-2.0, -1.0, 1.0) == -1.0, "mjolnir::compiletime::clamp");

namespace test
{
template<typename T>
decltype(boost::test_tools::tolerance(std::declval<T>())) tolerance();

template<>
decltype(boost::test_tools::tolerance(std::declval<float>()))
tolerance<float>()
{return boost::test_tools::tolerance(3.0f / static_cast<float>(std::pow(2, 12)));}

template<>
decltype(boost::test_tools::tolerance(std::declval<double>()))
tolerance<double>()
{return boost::test_tools::tolerance(2.0 / std::pow(2, 14));}
} // test

BOOST_AUTO_TEST_CASE_TEMPLATE(rsqrt, Real, test_targets)
{
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<Real> uni(0.0001, 100.0);

    for(std::size_t i=0; i < N; ++i)
    {
        const Real x     = uni(mt);
        const Real x_sq  = x * x;
        const Real x_inv = 1 / x;
        BOOST_TEST(x_inv == mjolnir::math::rsqrt(x_sq), test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(clamp, Real, test_targets)
{
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<Real> uni(-100.0, 100.0);
    const Real maximum =  10.0;
    const Real minimum = -10.0;

    for(std::size_t i=0; i < N; ++i)
    {
        const Real x = uni(mt);
        BOOST_TEST(mjolnir::math::clamp<Real>(x, minimum, maximum) ==
                   std::min(std::max(minimum, x), maximum),
                   test::tolerance<Real>());
    }
}

