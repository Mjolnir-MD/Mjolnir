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
