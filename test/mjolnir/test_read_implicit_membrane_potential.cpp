#define BOOST_TEST_MODULE "test_read_implicit_membrane_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_external_potential.hpp>

BOOST_AUTO_TEST_CASE(read_implicit_membrane_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_implicit_membrane.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"interaction",       toml::value("Distance")},
            {"potential",         toml::value("Excludedvolume_wallWall")},
            {"shape", toml::value(toml::table{
                    {"name",     toml::value("AxisAlignedPlane")},
                    {"axis",     toml::value("-X")},
                    {"position", toml::value(1.0)},
                    {"margin",   toml::value(0.5)}
            })},
            {"thickness",             toml::value(3.14)},
            {"interaction_magnitude", toml::value(6.28)},
            {"bend",                  toml::value(9.42)},
            {"parameters",            toml::value(toml::array{
                toml::table{{"index", 0}, {"hydrophobicity", 2.0}},
                toml::table{{"index", 1}, {"hydrophobicity", 2.0}}
            })}
        };
        const auto g = mjolnir::read_implicit_membrane_potential<real_type>(v);

        BOOST_TEST(g.hydrophobicities().size()     == 2);
        BOOST_TEST(g.hydrophobicities().at(0)      == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.hydrophobicities().at(1)      == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.half_thick()            == 3.14 * 0.5, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.interaction_magnitude() == 6.28, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.bend()                  == 9.42, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_implicit_membrane_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_implicit_membrane.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        const toml::value v = toml::table{
            {"interaction",       toml::value("Distance")},
            {"potential",         toml::value("Excludedvolume_wallWall")},
            {"shape", toml::value(toml::table{
                    {"name",     toml::value("AxisAlignedPlane")},
                    {"axis",     toml::value("-X")},
                    {"position", toml::value(1.0)},
                    {"margin",   toml::value(0.5)}
            })},
            {"thickness",             toml::value(3.14)},
            {"interaction_magnitude", toml::value(6.28)},
            {"bend",                  toml::value(9.42)},
            {"parameters",            toml::value(toml::array{
                toml::table{{"index", 0}, {"hydrophobicity", 2.0}},
                toml::table{{"index", 1}, {"hydrophobicity", 2.0}}
            })}
        };
        const auto g = mjolnir::read_implicit_membrane_potential<real_type>(v);

        BOOST_TEST(g.hydrophobicities().size()     == 2);
        BOOST_TEST(g.hydrophobicities().at(0)      == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.hydrophobicities().at(1)      == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.half_thick()             == 3.14f * 0.5f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.interaction_magnitude() == 6.28f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.bend()                  == 9.42f, boost::test_tools::tolerance(tol));
    }
}
