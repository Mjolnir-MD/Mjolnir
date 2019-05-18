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
        using namespace toml::literals;
        const auto v = u8R"(
            interaction    = "Distance"
            potential      = "ExcludedVolumeWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "-X"
            shape.position = 1.0
            shape.margin   = 0.5
            thickness      = 3.14
            interaction_magnitude = 6.28
            bend                  = 9.42
            parameters = [
                {index = 0, hydrophobicity = 2.0},
                {index = 1, hydrophobicity = 2.0},
            ]
        )"_toml;

        const auto g = mjolnir::read_implicit_membrane_potential<real_type>(v);

        BOOST_TEST(g.hydrophobicities().size()     == 2u);
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
        using namespace toml::literals;
        const auto v = u8R"(
            interaction    = "Distance"
            potential      = "ExcludedVolumeWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "-X"
            shape.position = 1.0
            shape.margin   = 0.5
            thickness      = 3.14
            interaction_magnitude = 6.28
            bend                  = 9.42
            parameters = [
                {index = 0, hydrophobicity = 2.0},
                {index = 1, hydrophobicity = 2.0},
            ]
        )"_toml;

        const auto g = mjolnir::read_implicit_membrane_potential<real_type>(v);

        BOOST_TEST(g.hydrophobicities().size()     == 2u);
        BOOST_TEST(g.hydrophobicities().at(0)      == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.hydrophobicities().at(1)      == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.half_thick()             == 3.14f * 0.5f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.interaction_magnitude() == 6.28f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.bend()                  == 9.42f, boost::test_tools::tolerance(tol));
    }
}
