#define BOOST_TEST_MODULE "test_read_periodic_gaussian_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_local_potential.hpp>

BOOST_AUTO_TEST_CASE(read_periodic_gaussian_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_periodic_gaussian.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",       toml::value(3.14)},
            {"sigma",   toml::value(.577)},
            {"v0",      toml::value(2.71)}
        };
        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }

    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",       toml::value(3.14)},
            {u8"σ",     toml::value(.577)},
            {"v0",      toml::value(2.71)}
        };
        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_periodic_gaussian_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_periodic_gaussian.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",       toml::value(3.14)},
            {"sigma",   toml::value(.577)},
            {"v0",      toml::value(2.71)}
        };
        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }

    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",       toml::value(3.14)},
            {u8"σ",     toml::value(.577)},
            {"v0",      toml::value(2.71)}
        };
        const auto g = mjolnir::read_periodic_gaussian_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.sigma() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }
}

// ---------------------------------------------------------------------------
// read_local_potential

BOOST_AUTO_TEST_CASE(read_local_potential_periodic_gaussian_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_periodic_gaussian.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"parameters", toml::value({toml::table{
                    {"indices", toml::value({1, 2})},
                    {"k",       toml::value(3.14)},
                    {"sigma",   toml::value(.577)},
                    {"v0",      toml::value(2.71)}
                }})}
        };
        const auto g = mjolnir::read_local_potential<2,
              mjolnir::PeriodicGaussianPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()     == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.sigma() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0()    == 2.71,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_local_potential_periodic_gaussian_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_periodic_gaussian.log");

    using real_type = float;
    constexpr real_type tol = 1e-4;
    {
        const toml::value v = toml::table{
            {"parameters", toml::value({toml::table{
                    {"indices", toml::value({1, 2})},
                    {"k",       toml::value(3.14)},
                    {"sigma",   toml::value(.577)},
                    {"v0",      toml::value(2.71)}
                }})}
        };
        const auto g = mjolnir::read_local_potential<2,
              mjolnir::PeriodicGaussianPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.sigma() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0()    == 2.71f,  boost::test_tools::tolerance(tol));
    }
}
