#define BOOST_TEST_MODULE "test_megingjord_simd_add"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include "utils.hpp"
#include <megingjord/simd/pack.hpp>

BOOST_AUTO_TEST_CASE(simd_add_tight_short)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    for(std::size_t count = 0; count < N; ++count)
    {
        const double x1 = uni(mt);
        const double x2 = uni(mt);
        const double x3 = uni(mt);
        const double x4 = uni(mt);

        const double y1 = uni(mt);
        const double y2 = uni(mt);
        const double y3 = uni(mt);
        const double y4 = uni(mt);

        const double z1 = x1 + y1;
        const double z2 = x2 + y2;
        const double z3 = x3 + y3;
        const double z4 = x4 + y4;

        megingjord::simd::packable_array<double, 4> xs{x1, x2, x3, x4};
        megingjord::simd::packable_array<double, 4> ys{y1, y2, y3, y4};
        megingjord::simd::packable_array<double, 4> zs = xs + ys;

        BOOST_CHECK_CLOSE_FRACTION(zs[0], z1, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(zs[1], z2, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(zs[2], z3, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(zs[3], z4, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(simd_add_loose_short)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    for(std::size_t count = 0; count < N; ++count)
    {
        const double x1 = uni(mt);
        const double x2 = uni(mt);
        const double x3 = uni(mt);

        const double y1 = uni(mt);
        const double y2 = uni(mt);
        const double y3 = uni(mt);

        const double z1 = x1 + y1;
        const double z2 = x2 + y2;
        const double z3 = x3 + y3;

        megingjord::simd::packable_array<double, 3> xs{x1, x2, x3};
        megingjord::simd::packable_array<double, 3> ys{y1, y2, y3};
        megingjord::simd::packable_array<double, 3> zs = xs + ys;

        BOOST_CHECK_CLOSE_FRACTION(zs[0], z1, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(zs[1], z2, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(zs[2], z3, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(simd_add_tight_long)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    for(std::size_t count = 0; count < N; ++count)
    {
        const auto rx = array_maker<double, tight>::invoke(mt, uni);
        const auto ry = array_maker<double, tight>::invoke(mt, uni);
        std::array<double, tight> rz;
        for(std::size_t i=0; i<tight; ++i)
            rz[i] = rx[i] + ry[i];

        megingjord::simd::packable_array<double, tight> xs;
        for(std::size_t i=0; i<tight; ++i) xs[i] = rx[i];

        megingjord::simd::packable_array<double, tight> ys;
        for(std::size_t i=0; i<tight; ++i) ys[i] = ry[i];

        megingjord::simd::packable_array<double, tight> zs = xs + ys;

        for(std::size_t i=0; i<tight; ++i)
            BOOST_CHECK_CLOSE_FRACTION(zs[i], rz[i], tolerance);
    }
}

BOOST_AUTO_TEST_CASE(simd_add_loose_long)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    for(std::size_t count = 0; count < N; ++count)
    {
        const auto rx = array_maker<double, loose>::invoke(mt, uni);
        const auto ry = array_maker<double, loose>::invoke(mt, uni);
        std::array<double, loose> rz;
        for(std::size_t i=0; i<loose; ++i)
            rz[i] = rx[i] + ry[i];

        megingjord::simd::packable_array<double, loose> xs;
        for(std::size_t i=0; i<loose; ++i) xs[i] = rx[i];

        megingjord::simd::packable_array<double, loose> ys;
        for(std::size_t i=0; i<loose; ++i) ys[i] = ry[i];

        megingjord::simd::packable_array<double, loose> zs = xs + ys;

        for(std::size_t i=0; i<loose; ++i)
            BOOST_CHECK_CLOSE_FRACTION(zs[i], rz[i], tolerance);
    }
}



