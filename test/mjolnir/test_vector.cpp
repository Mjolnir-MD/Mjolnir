#define BOOST_TEST_MODULE "test_vector"

#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include <mjolnir/math/Vector.hpp>

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
{return boost::test_tools::tolerance(3.0f / static_cast<float>(std::pow(2, 12)));}

template<>
decltype(boost::test_tools::tolerance(std::declval<double>()))
tolerance<double>()
{return boost::test_tools::tolerance(2.0 / std::pow(2, 14));}
} // test

BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_add, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-100.0, 100.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const math::Vector<Real, 3> lhs(uni(mt), uni(mt), uni(mt)),
                              rhs(uni(mt), uni(mt), uni(mt));
        const auto add = lhs + rhs;
        BOOST_TEST(add[0] == lhs[0] + rhs[0], test::tolerance<Real>());
        BOOST_TEST(add[1] == lhs[1] + rhs[1], test::tolerance<Real>());
        BOOST_TEST(add[2] == lhs[2] + rhs[2], test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_sub, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-100.0, 100.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const math::Vector<Real, 3> lhs(uni(mt), uni(mt), uni(mt)),
                              rhs(uni(mt), uni(mt), uni(mt));
        const auto sub = lhs - rhs;
        BOOST_TEST(sub[0] == lhs[0] - rhs[0], test::tolerance<Real>());
        BOOST_TEST(sub[1] == lhs[1] - rhs[1], test::tolerance<Real>());
        BOOST_TEST(sub[2] == lhs[2] - rhs[2], test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_mul, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-100.0, 100.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const Real            lhs(uni(mt));
        const math::Vector<Real, 3> rhs(uni(mt), uni(mt), uni(mt));
        const auto mul = lhs * rhs;
        BOOST_TEST(mul[0] == lhs * rhs[0], test::tolerance<Real>());
        BOOST_TEST(mul[1] == lhs * rhs[1], test::tolerance<Real>());
        BOOST_TEST(mul[2] == lhs * rhs[2], test::tolerance<Real>());
    }

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const math::Vector<Real, 3> lhs(uni(mt), uni(mt), uni(mt));
        const Real            rhs(uni(mt));
        const auto mul = lhs * rhs;
        BOOST_TEST(mul[0] == lhs[0] * rhs, test::tolerance<Real>());
        BOOST_TEST(mul[1] == lhs[1] * rhs, test::tolerance<Real>());
        BOOST_TEST(mul[2] == lhs[2] * rhs, test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_div, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-100.0, 100.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const math::Vector<Real, 3> lhs(uni(mt), uni(mt), uni(mt));
        const Real            rhs(uni(mt));
        const auto div = lhs / rhs;
        BOOST_TEST(div[0] == lhs[0] / rhs, test::tolerance<Real>());
        BOOST_TEST(div[1] == lhs[1] / rhs, test::tolerance<Real>());
        BOOST_TEST(div[2] == lhs[2] / rhs, test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_dot, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const math::Vector<Real, 3> lhs(uni(mt), uni(mt), uni(mt));
        const math::Vector<Real, 3> rhs(uni(mt), uni(mt), uni(mt));
        const Real dot = math::dot_product(lhs, rhs);
        BOOST_TEST(dot == lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2],
                   test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_len, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const math::Vector<Real, 3> lhs(uni(mt), uni(mt), uni(mt));
        const Real dot   = math::dot_product(lhs, lhs);
        const Real lensq = math::length_sq(lhs);
        const Real len   = math::length(lhs);
        BOOST_TEST(lensq == dot, test::tolerance<Real>());

        BOOST_TEST(lensq == dot,            test::tolerance<Real>());
        BOOST_TEST(len   == std::sqrt(dot), test::tolerance<Real>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_vector_cross_product, Real, test_targets)
{
    using namespace mjolnir;

    std::mt19937 mt(seed);
    std::uniform_real_distribution<Real> uni(-1.0, 1.0);

    for(std::size_t test_times=0; test_times<N; ++test_times)
    {
        const math::Vector<Real, 3> lhs(uni(mt), uni(mt), uni(mt));
        const math::Vector<Real, 3> rhs(uni(mt), uni(mt), uni(mt));
        const math::Vector<Real, 3> cross = math::cross_product(lhs, rhs);
        const Real dotl = math::dot_product(cross, lhs);
        const Real dotr = math::dot_product(cross, rhs);

        BOOST_TEST(dotl == static_cast<Real>(0.0), test::tolerance<Real>());
        BOOST_TEST(dotr == static_cast<Real>(0.0), test::tolerance<Real>());
        const Real lenc = math::length(cross);

        const Real lenl = math::length(lhs);
        const Real lenr = math::length(rhs);
        const Real dot  = math::dot_product(lhs, rhs);
        const Real cost = dot / (lenl * lenr);
        const Real sint = std::sqrt(1.0 - cost * cost);

        BOOST_TEST(lenc == lenl * lenr * sint, test::tolerance<Real>());
    }
}
