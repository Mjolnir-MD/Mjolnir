#define BOOST_TEST_MODULE "test_read_integrator"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_integrator.hpp>

BOOST_AUTO_TEST_CASE(read_velocity_verlet_integrator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_integrator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{
            {"delta_t",       toml::value(0.1)},
        };

        const auto integr = mjolnir::read_velocity_verlet_integrator<traits_type>(v);
        BOOST_TEST(integr.delta_t() == 0.1, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_underdamped_langevin_integrator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_integrator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{
            {"delta_t",       toml::value(0.1)},
            {"seed",          toml::value(1234)},
            {"parameters",    toml::value(toml::array{
                toml::table{{"index", toml::value(0)}, {"gamma", toml::value(0.1)}},
                toml::table{{"index", toml::value(1)}, {u8"γ",   toml::value(0.2)}},
            })},
        };

        const auto integr = mjolnir::read_underdamped_langevin_integrator<traits_type>(v);
        BOOST_TEST(integr.delta_t() == 0.1, boost::test_tools::tolerance(tol));
        BOOST_TEST(integr.parameters().size() == 2);
        BOOST_TEST(integr.parameters().at(0) == 0.1, boost::test_tools::tolerance(tol));
        BOOST_TEST(integr.parameters().at(1) == 0.2, boost::test_tools::tolerance(tol));
    }
}
