#define BOOST_TEST_MODULE "test_read_excluded_volume_wall_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_excluded_volume_wall_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_excluded_volume_wall.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction    = "Distance"
            potential      = "ExcludedVolumeWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "-X"
            shape.position = 1.0
            shape.margin   = 0.5
            epsilon        = 3.14
            parameters     = [
                {index = 0, radius = 2.0},
                {index = 1, radius = 3.0},
                {index = 5, radius = 4.0},
            ]
        )"_toml;

        const auto g = mjolnir::read_excluded_volume_wall_potential<real_type>(v);

        BOOST_TEST(g.participants().size() == 3u);
        BOOST_TEST(g.participants().at(0) == 0u);
        BOOST_TEST(g.participants().at(1) == 1u);
        BOOST_TEST(g.participants().at(2) == 5u);
        BOOST_TEST(g.parameters().at(0)  == real_type(2.0 ), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(1)  == real_type(3.0 ), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(5)  == real_type(4.0 ), tolerance<real_type>());
        BOOST_TEST(g.epsilon()           == real_type(3.14), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_excluded_volume_wall_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_excluded_volume_wall.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction    = "Distance"
            potential      = "ExcludedVolumeWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "-X"
            shape.position = 1.0
            shape.margin   = 0.5
            epsilon        = 3.14
            env.large = 10.0
            env.small = 0.01
            parameters     = [
                {index = 0, radius = "large"},
                {index = 1, radius = "small"},
                {index = 5, radius = "large"},
            ]
        )"_toml;

        const auto g = mjolnir::read_excluded_volume_wall_potential<real_type>(v);

        BOOST_TEST(g.participants().size() == 3u);
        BOOST_TEST(g.participants().at(0) == 0u);
        BOOST_TEST(g.participants().at(1) == 1u);
        BOOST_TEST(g.participants().at(2) == 5u);

        BOOST_TEST(g.parameters().at(0)  == real_type(10.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(1)  == real_type(0.01), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(5)  == real_type(10.0), tolerance<real_type>());
        BOOST_TEST(g.epsilon()           == real_type(3.14), tolerance<real_type>());
    }
}
