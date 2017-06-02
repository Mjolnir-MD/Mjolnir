#define BOOST_TEST_MODULE "test_fast_inv_sqrt"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <cmath>

#include <random>
#include <chrono>
#include <array>
#include <iostream>
constexpr static unsigned int seed = 32479327;
constexpr static std::size_t N = 10000;

BOOST_AUTO_TEST_CASE(fast_inv_sqrt_double)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<double> uni(0.0, 100.0);
    for(std::size_t i=0; i < N; ++i)
    {
        const double x = uni(mt);
        const double x2 = x * x;
        const double x_inv = 1. / x;
        BOOST_CHECK_CLOSE_FRACTION(x_inv, mjolnir::fast_inv_sqrt(x2), 1e-10);
    }
}

BOOST_AUTO_TEST_CASE(fast_inv_sqrt_float)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<float> uni(0.0, 100.0);
    for(std::size_t i=0; i < N; ++i)
    {
        const float x = uni(mt);
        const float x2 = x * x;
        const float x_inv = 1. / x;
        BOOST_CHECK_CLOSE_FRACTION(x_inv, mjolnir::fast_inv_sqrt(x2), 1e-5);
    }
}

// BOOST_AUTO_TEST_CASE(speed_check_double)
// {
//     std::mt19937 mt(seed);
//     std::uniform_real_distribution<double> uni(0.0, 100.0);
//     std::array<double, N> rnds;
//     for(std::size_t i=0; i < N; ++i)
//         rnds[i] = uni(mt);
//
//     const auto normal_start = std::chrono::system_clock::now();
//     std::array<double, N> tmp1;
//     for(std::size_t i=0; i<N; ++i)
//         tmp1[i] = 1. / std::sqrt(rnds[i]);
//     const auto normal_end = std::chrono::system_clock::now();
//     const auto normal_duration = normal_end - normal_start;
//
//     const auto faster_start = std::chrono::system_clock::now();
//     std::array<double, N> tmp2;
//     for(std::size_t i=0; i<N; ++i)
//         tmp2[i] = mjolnir::fast_inv_sqrt(rnds[i]);
//     const auto faster_end = std::chrono::system_clock::now();
//     const auto faster_duration = faster_end - faster_start;
//
//     BOOST_CHECK(faster_duration < normal_duration);
// }
//
//
// BOOST_AUTO_TEST_CASE(speed_check_float)
// {
//     std::mt19937 mt(seed);
//     std::uniform_real_distribution<float> uni(0.0, 100.0);
//     std::array<float, N> rnds;
//     for(std::size_t i=0; i < N; ++i)
//         rnds[i] = uni(mt);
//
//     const auto normal_start = std::chrono::system_clock::now();
//     std::array<float, N> tmp1;
//     for(std::size_t i=0; i<N; ++i)
//         tmp1[i] = 1. / std::sqrt(rnds[i]);
//     const auto normal_end = std::chrono::system_clock::now();
//     const auto normal_duration = normal_end - normal_start;
//
//     const auto faster_start = std::chrono::system_clock::now();
//     std::array<float, N> tmp2;
//     for(std::size_t i=0; i<N; ++i)
//         tmp2[i] = mjolnir::fast_inv_sqrt(rnds[i]);
//     const auto faster_end = std::chrono::system_clock::now();
//     const auto faster_duration = faster_end - faster_start;
//
//     BOOST_CHECK(faster_duration < normal_duration);
//
// }
