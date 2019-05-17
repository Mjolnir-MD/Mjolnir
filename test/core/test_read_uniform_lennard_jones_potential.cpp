#define BOOST_TEST_MODULE "test_read_uniform_lennard_jones_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_global_potential.hpp>

BOOST_AUTO_TEST_CASE(read_uniform_lennard_jones_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_lennard_jones.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            sigma   = 2.0
            epsilon = 1.5
        )"_toml;

        const auto g = mjolnir::read_uniform_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.sigma()   == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon() == 1.5,  boost::test_tools::tolerance(tol));
    }
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            "σ" = 2.0
            "ε" = 1.5
        )"_toml;

        const auto g = mjolnir::read_uniform_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.sigma()   == 2.0,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon() == 1.5,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_uniform_lennard_jones_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_lennard_jones.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            sigma   = 2.0
            epsilon = 1.5
        )"_toml;

        const auto g = mjolnir::read_uniform_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.sigma()   == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon() == 1.5f,  boost::test_tools::tolerance(tol));
    }
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            "σ" = 2.0
            "ε" = 1.5
        )"_toml;

        const auto g = mjolnir::read_uniform_lennard_jones_potential<real_type>(v);

        const auto ignore_within = g.ignore_within();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(g.ignore_within().size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(g.sigma()   == 2.0f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.epsilon() == 1.5f,  boost::test_tools::tolerance(tol));
    }
}
