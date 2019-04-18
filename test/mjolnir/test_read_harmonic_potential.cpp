#define BOOST_TEST_MODULE "test_read_harmonic_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_local_potential.hpp>

BOOST_AUTO_TEST_CASE(read_harmonic_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
            v0      = 2.71
        )"_toml;

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
        using namespace toml::literals;
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
            v0      = 2.71
        )"_toml;

        const auto g = mjolnir::read_harmonic_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_local_potential_harmonic_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            parameters = [
                {indices = [1,2], k = 3.14, v0 = 2.71}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
            mjolnir::HarmonicPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()  == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0() == 2.71,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_local_potential_harmonic_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        using namespace toml::literals;
        const auto v = u8R"(
            parameters = [
                {indices = [1,2], k = 3.14, v0 = 2.71}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
            mjolnir::HarmonicPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()  == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0() == 2.71f,  boost::test_tools::tolerance(tol));
    }
}
