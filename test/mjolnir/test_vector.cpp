#define BOOST_TEST_MODULE "test_vector"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/math/Vector.hpp>

#include <random>
using namespace mjolnir;
typedef double real_type;
typedef mjolnir::Vector<double, 3> Vector3d;
constexpr static unsigned int seed = 32479327;
constexpr static std::size_t N = 10000;
constexpr static real_type tolerance = 1e-8;

BOOST_AUTO_TEST_CASE(test_vector_add)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Vector3d lhs(uni(mt), uni(mt), uni(mt));
        const Vector3d rhs(uni(mt), uni(mt), uni(mt));
        const Vector3d sum = lhs + rhs;
        BOOST_CHECK_CLOSE_FRACTION(sum[0], lhs[0] + rhs[0], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(sum[1], lhs[1] + rhs[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(sum[2], lhs[2] + rhs[2], tolerance);
    }
}

BOOST_AUTO_TEST_CASE(test_vector_subtract)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Vector3d lhs(uni(mt), uni(mt), uni(mt));
        const Vector3d rhs(uni(mt), uni(mt), uni(mt));
        const Vector3d sub = lhs - rhs;
        BOOST_CHECK_CLOSE_FRACTION(sub[0], lhs[0] - rhs[0], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(sub[1], lhs[1] - rhs[1], tolerance);
        BOOST_CHECK_CLOSE_FRACTION(sub[2], lhs[2] - rhs[2], tolerance);
    }
}

BOOST_AUTO_TEST_CASE(test_vector_scalar_multiple)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Vector3d lhs(uni(mt), uni(mt), uni(mt));
        const real_type rhs = uni(mt);
        const Vector3d mul = lhs * rhs;
        BOOST_CHECK_CLOSE_FRACTION(mul[0], lhs[0] * rhs, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(mul[1], lhs[1] * rhs, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(mul[2], lhs[2] * rhs, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(test_vector_scalar_division)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Vector3d lhs(uni(mt), uni(mt), uni(mt));
        const real_type rhs = uni(mt);
        const Vector3d div = lhs / rhs;
        BOOST_CHECK_CLOSE_FRACTION(div[0], lhs[0] / rhs, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(div[1], lhs[1] / rhs, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(div[2], lhs[2] / rhs, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(test_vector_dot_product)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Vector3d lhs(uni(mt), uni(mt), uni(mt));
        const Vector3d rhs(uni(mt), uni(mt), uni(mt));
        const real_type dot = dot_product(lhs, rhs);
        BOOST_CHECK_CLOSE_FRACTION(dot, lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2], tolerance);
    }
}

BOOST_AUTO_TEST_CASE(test_vector_length)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Vector3d lhs(uni(mt), uni(mt), uni(mt));
        const real_type dot = dot_product(lhs, lhs);
        const real_type lensq = length_sq(lhs);
        BOOST_CHECK_CLOSE_FRACTION(lensq, dot, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(test_vector_cross_product)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Vector3d lhs(uni(mt), uni(mt), uni(mt));
        const Vector3d rhs(uni(mt), uni(mt), uni(mt));
        const Vector3d cross = cross_product(lhs, rhs);
        const real_type dotl = dot_product(cross, lhs);
        const real_type dotr = dot_product(cross, rhs);
        BOOST_CHECK_SMALL(dotl, tolerance);
        BOOST_CHECK_SMALL(dotr, tolerance);
        const real_type lenc = length(cross);

        const real_type lenl = length(lhs);
        const real_type lenr = length(rhs);
        const real_type dot  = dot_product(lhs, rhs);
        const real_type cost = dot / (lenl * lenr);
        const real_type sint = std::sqrt(1. - cost * cost);

        BOOST_CHECK_CLOSE_FRACTION(lenc, lenl * lenr * sint, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(test_vector_scalar_triple_product)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Vector3d lhs(uni(mt), uni(mt), uni(mt));
        const Vector3d mid(uni(mt), uni(mt), uni(mt));
        const Vector3d rhs(uni(mt), uni(mt), uni(mt));
        const real_type tri = scalar_triple_product(lhs, mid, rhs);
        BOOST_CHECK_CLOSE_FRACTION(tri, dot_product(cross_product(lhs, mid), rhs), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(tri, dot_product(lhs, cross_product(mid, rhs)), tolerance);
    }
}

BOOST_AUTO_TEST_CASE(test_vector_rotation)
{
    std::mt19937 mt(seed);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    const Vector3d unit_x(1., 0., 0.);
    const Vector3d unit_z(0., 0., 1.);
    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const real_type theta = uni(mt) * M_PI;
        const Vector3d rot = rotate(theta, unit_z, unit_x);
        BOOST_CHECK_CLOSE_FRACTION(length(rot), 1., tolerance);
        BOOST_CHECK_CLOSE_FRACTION(rot[0], std::cos(theta), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(rot[1], std::sin(theta), tolerance);
    }
}


