#define BOOST_TEST_MODULE "test_read_harmonic_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/input/read_potential.hpp>

BOOST_AUTO_TEST_CASE(read_harmonic_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",       toml::value(3.14)},
            {"v0",      toml::value(2.71)},
        };
        const auto g = mjolnir::read_harmonic_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_harmonic_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",       toml::value(3.14)},
            {"v0",      toml::value(2.71)},
        };
        const auto g = mjolnir::read_harmonic_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }
}
