#define BOOST_TEST_MODULE "test_read_excluded_volume_wall_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/input/read_external_potential.hpp>

BOOST_AUTO_TEST_CASE(read_excluded_volume_wall_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_excluded_volume_wall.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"interaction",       toml::value("Distance")},
            {"potential",         toml::value("ExcludedvolumeWall")},
            {"shape", toml::value(toml::table{
                    {"name",     toml::value("AxisAlignedPlane")},
                    {"axis",     toml::value("-X")},
                    {"position", toml::value(1.0)},
                    {"margin",   toml::value(0.5)}
            })},
            {"epsilon",           toml::value(3.14)},
            {"parameters",        toml::value(toml::array{
                toml::table{{"index", 0}, {"radius", 2.0}},
                toml::table{{"index", 1}, {"radius", 2.0}}
            })}
        };
        const auto g = mjolnir::read_excluded_volume_wall_potential<real_type>(v);

        BOOST_TEST(g.parameters().size() == 2);
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
        const toml::value v = toml::table{
            {"interaction",       toml::value("Distance")},
            {"potential",         toml::value("ExcludedVolumeWall")},
            {"shape", toml::value(toml::table{
                    {"name",     toml::value("AxisAlignedPlane")},
                    {"axis",     toml::value("-X")},
                    {"position", toml::value(1.0)},
                    {"margin",   toml::value(0.5)}
            })},
            {"epsilon",           toml::value(3.14)},
            {"parameters",        toml::value(toml::array{
                toml::table{{"index", 0}, {"radius", 2.0}},
                toml::table{{"index", 1}, {"radius", 2.0}}
            })}
        };
        const auto g = mjolnir::read_excluded_volume_wall_potential<real_type>(v);

        BOOST_TEST(g.parameters().size() == 2);
        BOOST_TEST(g.parameters().at(0)  == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1)  == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon()           == 3.14f, boost::test_tools::tolerance(tol));
    }
}
