#define BOOST_TEST_MODULE "test_read_excluded_volume_wall_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_external_potential.hpp>

BOOST_AUTO_TEST_CASE(read_excluded_volume_wall_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_excluded_volume_wall.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
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
                {index = 1, radius = 2.0},
            ]
        )"_toml;

        const auto g = mjolnir::read_excluded_volume_wall_potential<real_type>(v);

        BOOST_TEST(g.parameters().size() == 2u);
        BOOST_TEST(g.parameters().at(0)  == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1)  == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon()           == 3.14, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_excluded_volume_wall_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_excluded_volume_wall.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction             = "Distance"
            potential               = "ExcludedVolumeWall"
            shape.name              = "AxisAlignedPlane"
            shape.axis              = "-X"
            shape.position          = 1.0
            shape.margin            = 0.5
            epsilon                 = 3.14
            parameters  = [
                {index = 0, radius = 2.0},
                {index = 1, radius = 2.0},
            ]
        )"_toml;

        const auto g = mjolnir::read_excluded_volume_wall_potential<real_type>(v);

        BOOST_TEST(g.parameters().size() == 2u);
        BOOST_TEST(g.parameters().at(0)  == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1)  == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon()           == 3.14f, boost::test_tools::tolerance(tol));
    }
}
