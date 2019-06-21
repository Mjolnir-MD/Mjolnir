#define BOOST_TEST_MODULE "test_read_lennard_jones_wall_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_external_potential.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_lennard_jones_wall_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_lennard_jones_wall.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction    = "Distance"
            potential      = "LennardJonesWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "+X"
            shape.position = 1.0
            shape.margin   = 0.5
            parameters     = [
                {index = 0, sigma = 2.0, epsilon = 1.5},
                {index = 1, "σ"   = 2.0, epsilon = 1.5},
                {index = 2, sigma = 2.0, "ε"     = 1.5},
                {index = 3, "σ"   = 2.0, "ε"     = 1.5},
            ]
        )"_toml;

        const auto g = mjolnir::read_lennard_jones_wall_potential<real_type>(v);

        BOOST_TEST(g.parameters().size() == 4u);
        BOOST_TEST(g.parameters().at(0).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(1).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(2).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(3).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(0).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(1).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(2).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(3).second == real_type(1.5), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_lennard_jones_wall_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_lennard_jones_wall.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction    = "Distance"
            potential      = "LennardJonesWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "+X"
            shape.position = 1.0
            shape.margin   = 0.5
            env.large  = 10.0
            env.small  = 0.01
            env.strong = 100.0
            env.weak   = 0.001
            parameters     = [
                {index = 0, sigma = "large", epsilon = "strong"},
                {index = 1, "σ"   = "large", epsilon = "weak"},
                {index = 2, sigma = "small", "ε"     = "strong"},
                {index = 3, "σ"   = "small", "ε"     = "weak"},
            ]
        )"_toml;

        const auto g = mjolnir::read_lennard_jones_wall_potential<real_type>(v);

        BOOST_TEST(g.parameters().size() == 4u);
        BOOST_TEST(g.parameters().at(0).first  == real_type(10.0),  tolerance<real_type>());
        BOOST_TEST(g.parameters().at(1).first  == real_type(10.0),  tolerance<real_type>());
        BOOST_TEST(g.parameters().at(2).first  == real_type(0.01),  tolerance<real_type>());
        BOOST_TEST(g.parameters().at(3).first  == real_type(0.01),  tolerance<real_type>());
        BOOST_TEST(g.parameters().at(0).second == real_type(100.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(1).second == real_type(0.001), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(2).second == real_type(100.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(3).second == real_type(0.001), tolerance<real_type>());
    }
}
