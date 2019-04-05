#define BOOST_TEST_MODULE "test_read_lennard_jones_wall_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_external_potential.hpp>

BOOST_AUTO_TEST_CASE(read_lennard_jones_wall_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_lennard_jones_wall.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"interaction",       toml::value("Distance")},
            {"potential",         toml::value("LennardJonesWall")},
            {"shape", toml::value(toml::table{
                    {"name",     toml::value("AxisAlignedPlane")},
                    {"axis",     toml::value("-X")},
                    {"position", toml::value(1.0)},
                    {"margin",   toml::value(0.5)}
            })},
            {"parameters",        toml::value(toml::array{
                toml::table{{"index", 0}, {"sigma", 2.0}, {"epsilon", 1.5}},
                toml::table{{"index", 1}, {u8"σ",   2.0}, {u8"ε",     1.5}},
                toml::table{{"index", 2}, {u8"σ",   2.0}, {"epsilon", 1.5}},
                toml::table{{"index", 3}, {"sigma", 2.0}, {u8"ε",     1.5}}
            })}
        };
        const auto g = mjolnir::read_lennard_jones_wall_potential<real_type>(v);

        BOOST_TEST(g.parameters().size() == 4);
        BOOST_TEST(g.parameters().at(0).first  == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1).first  == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(2).first  == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(3).first  == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(0).second == 1.5,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1).second == 1.5,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(2).second == 1.5,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(3).second == 1.5,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_lennard_jones_wall_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_lennard_jones_wall.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        const toml::value v = toml::table{
            {"interaction",       toml::value("Distance")},
            {"potential",         toml::value("LennardJonesWall")},
            {"shape", toml::value(toml::table{
                    {"name",     toml::value("AxisAlignedPlane")},
                    {"axis",     toml::value("-X")},
                    {"position", toml::value(1.0)},
                    {"margin",   toml::value(0.5)}
            })},
            {"parameters",        toml::value(toml::array{
                toml::table{{"index", 0}, {"sigma", 2.0}, {"epsilon", 1.5}},
                toml::table{{"index", 1}, {u8"σ",   2.0}, {u8"ε",     1.5}},
                toml::table{{"index", 2}, {u8"σ",   2.0}, {"epsilon", 1.5}},
                toml::table{{"index", 3}, {"sigma", 2.0}, {u8"ε",     1.5}}
            })}
        };
        const auto g = mjolnir::read_lennard_jones_wall_potential<real_type>(v);

        BOOST_TEST(g.parameters().size() == 4);
        BOOST_TEST(g.parameters().at(0).first  == 2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1).first  == 2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(2).first  == 2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(3).first  == 2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(0).second == 1.5f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(1).second == 1.5f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(2).second == 1.5f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.parameters().at(3).second == 1.5f, boost::test_tools::tolerance(tol));
    }
}
