#define BOOST_TEST_MODULE "test_matrix"

#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include <mjolnir/math/Matrix.hpp>
#include <random>
#include <cstdint>

constexpr std::uint32_t seed = 123456789;
constexpr std::size_t   N    = 10000;
typedef boost::mpl::list<double, float> test_targets;

namespace test
{
template<typename T>
decltype(boost::test_tools::tolerance(std::declval<T>())) tolerance();

template<>
decltype(boost::test_tools::tolerance(std::declval<float>()))
tolerance<float>()
{return boost::test_tools::tolerance(3.0f / static_cast<float>(std::pow(2, 8)));}

template<>
decltype(boost::test_tools::tolerance(std::declval<double>()))
tolerance<double>()
{return boost::test_tools::tolerance(2.0 / std::pow(2, 14));}
} // test

BOOST_AUTO_TEST_CASE_TEMPLATE(add_matrix_3x3, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Matrix<Real, 3, 3> lhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            ), rhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            );
        const auto add = lhs + rhs;

        BOOST_TEST(add(0, 0) == lhs(0, 0) + rhs(0, 0), test::tolerance<Real>());
        BOOST_TEST(add(0, 1) == lhs(0, 1) + rhs(0, 1), test::tolerance<Real>());
        BOOST_TEST(add(0, 2) == lhs(0, 2) + rhs(0, 2), test::tolerance<Real>());
        BOOST_TEST(add(1, 0) == lhs(1, 0) + rhs(1, 0), test::tolerance<Real>());
        BOOST_TEST(add(1, 1) == lhs(1, 1) + rhs(1, 1), test::tolerance<Real>());
        BOOST_TEST(add(1, 2) == lhs(1, 2) + rhs(1, 2), test::tolerance<Real>());
        BOOST_TEST(add(2, 0) == lhs(2, 0) + rhs(2, 0), test::tolerance<Real>());
        BOOST_TEST(add(2, 1) == lhs(2, 1) + rhs(2, 1), test::tolerance<Real>());
        BOOST_TEST(add(2, 2) == lhs(2, 2) + rhs(2, 2), test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sub_matrix_3x3, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Matrix<Real, 3, 3> lhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            ), rhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            );
        const auto sub = lhs - rhs;

        BOOST_TEST(sub(0, 0) == lhs(0, 0) - rhs(0, 0), test::tolerance<Real>());
        BOOST_TEST(sub(0, 1) == lhs(0, 1) - rhs(0, 1), test::tolerance<Real>());
        BOOST_TEST(sub(0, 2) == lhs(0, 2) - rhs(0, 2), test::tolerance<Real>());
        BOOST_TEST(sub(1, 0) == lhs(1, 0) - rhs(1, 0), test::tolerance<Real>());
        BOOST_TEST(sub(1, 1) == lhs(1, 1) - rhs(1, 1), test::tolerance<Real>());
        BOOST_TEST(sub(1, 2) == lhs(1, 2) - rhs(1, 2), test::tolerance<Real>());
        BOOST_TEST(sub(2, 0) == lhs(2, 0) - rhs(2, 0), test::tolerance<Real>());
        BOOST_TEST(sub(2, 1) == lhs(2, 1) - rhs(2, 1), test::tolerance<Real>());
        BOOST_TEST(sub(2, 2) == lhs(2, 2) - rhs(2, 2), test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(scalar_mul_matrix_3x3, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Real lhs = uni(mt);
        const Matrix<Real, 3, 3> rhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            );
        const auto mul = lhs * rhs;

        BOOST_TEST(mul(0, 0) == lhs * rhs(0, 0), test::tolerance<Real>());
        BOOST_TEST(mul(0, 1) == lhs * rhs(0, 1), test::tolerance<Real>());
        BOOST_TEST(mul(0, 2) == lhs * rhs(0, 2), test::tolerance<Real>());
        BOOST_TEST(mul(1, 0) == lhs * rhs(1, 0), test::tolerance<Real>());
        BOOST_TEST(mul(1, 1) == lhs * rhs(1, 1), test::tolerance<Real>());
        BOOST_TEST(mul(1, 2) == lhs * rhs(1, 2), test::tolerance<Real>());
        BOOST_TEST(mul(2, 0) == lhs * rhs(2, 0), test::tolerance<Real>());
        BOOST_TEST(mul(2, 1) == lhs * rhs(2, 1), test::tolerance<Real>());
        BOOST_TEST(mul(2, 2) == lhs * rhs(2, 2), test::tolerance<Real>());
    }

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Matrix<Real, 3, 3> lhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            );
        const Real rhs = uni(mt);
        const auto mul = lhs * rhs;

        BOOST_TEST(mul(0, 0) == lhs(0, 0) * rhs, test::tolerance<Real>());
        BOOST_TEST(mul(0, 1) == lhs(0, 1) * rhs, test::tolerance<Real>());
        BOOST_TEST(mul(0, 2) == lhs(0, 2) * rhs, test::tolerance<Real>());
        BOOST_TEST(mul(1, 0) == lhs(1, 0) * rhs, test::tolerance<Real>());
        BOOST_TEST(mul(1, 1) == lhs(1, 1) * rhs, test::tolerance<Real>());
        BOOST_TEST(mul(1, 2) == lhs(1, 2) * rhs, test::tolerance<Real>());
        BOOST_TEST(mul(2, 0) == lhs(2, 0) * rhs, test::tolerance<Real>());
        BOOST_TEST(mul(2, 1) == lhs(2, 1) * rhs, test::tolerance<Real>());
        BOOST_TEST(mul(2, 2) == lhs(2, 2) * rhs, test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(scalar_div_matrix_3x3, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Matrix<Real, 3, 3> lhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            );
        const Real rhs = uni(mt);
        const auto mul = lhs / rhs;

        BOOST_TEST(mul(0, 0) == lhs(0, 0) / rhs, test::tolerance<Real>());
        BOOST_TEST(mul(0, 1) == lhs(0, 1) / rhs, test::tolerance<Real>());
        BOOST_TEST(mul(0, 2) == lhs(0, 2) / rhs, test::tolerance<Real>());
        BOOST_TEST(mul(1, 0) == lhs(1, 0) / rhs, test::tolerance<Real>());
        BOOST_TEST(mul(1, 1) == lhs(1, 1) / rhs, test::tolerance<Real>());
        BOOST_TEST(mul(1, 2) == lhs(1, 2) / rhs, test::tolerance<Real>());
        BOOST_TEST(mul(2, 0) == lhs(2, 0) / rhs, test::tolerance<Real>());
        BOOST_TEST(mul(2, 1) == lhs(2, 1) / rhs, test::tolerance<Real>());
        BOOST_TEST(mul(2, 2) == lhs(2, 2) / rhs, test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(mul_matrix_3x3, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Matrix<Real, 3, 3> lhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            ), rhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            );
        const auto mul = lhs * rhs;

        for(std::size_t i=0; i<3; ++i)
        {
            for(std::size_t j=0; j<3; ++j)
            {
                Real sum = 0.;
                for(std::size_t k=0; k<3; ++k)
                {
                    sum += lhs(i, k) * rhs(k, j);
                }
                BOOST_TEST(mul(i, j) == sum, test::tolerance<Real>());
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(matrix_3x3_inverse, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Matrix<Real, 3, 3> lhs(
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt),
                uni(mt), uni(mt), uni(mt)
            );
        if(determinant(lhs) == 0.0)
        {
            continue;
        }

        const auto inv = inverse(lhs);
        const auto unit1 = lhs * inv;
        const auto unit2 = inv * lhs;

        BOOST_TEST(unit1(0, 0) == static_cast<Real>(1.0), test::tolerance<Real>());
        BOOST_TEST(unit1(1, 1) == static_cast<Real>(1.0), test::tolerance<Real>());
        BOOST_TEST(unit1(2, 2) == static_cast<Real>(1.0), test::tolerance<Real>());
        BOOST_TEST(unit1(0, 1) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit1(0, 2) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit1(1, 0) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit1(1, 2) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit1(2, 0) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit1(2, 1) == static_cast<Real>(0.0), test::tolerance<Real>());

        BOOST_TEST(unit2(0, 0) == static_cast<Real>(1.0), test::tolerance<Real>());
        BOOST_TEST(unit2(1, 1) == static_cast<Real>(1.0), test::tolerance<Real>());
        BOOST_TEST(unit2(2, 2) == static_cast<Real>(1.0), test::tolerance<Real>());
        BOOST_TEST(unit2(0, 1) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit2(0, 2) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit2(1, 0) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit2(1, 2) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit2(2, 0) == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(unit2(2, 1) == static_cast<Real>(0.0), test::tolerance<Real>());
    }
}

