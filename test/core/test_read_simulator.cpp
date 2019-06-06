#define BOOST_TEST_MODULE "test_read_simulator"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <test/util/make_empty_input.hpp>
#include <mjolnir/input/read_simulator.hpp>

BOOST_AUTO_TEST_CASE(read_newtonian_molecular_dynamics_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_simulator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    auto root = mjolnir::test::make_empty_input();
    {
        using namespace toml::literals;
        const auto v = u8R"(
            type            = "MolecularDynamics"
            integrator.type = "VelocityVerlet"
            precision       = "double"
            boundary_type   = "Unlimited"
            delta_t         = 0.1
            total_step      = 100
            save_step       = 10
        )"_toml;

        root["simulator"] = v;
        const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto mdsim = dynamic_cast<mjolnir::MolecularDynamicsSimulator<
            traits_type, mjolnir::VelocityVerletIntegrator<traits_type>>*>(sim.get());
        BOOST_TEST(static_cast<bool>(mdsim));

        sim->initialize();
        for(std::size_t i=0; i<99; ++i)
        {
            BOOST_TEST(mdsim->time() == i * 0.1, boost::test_tools::tolerance(tol));
            BOOST_TEST(sim->step()); // check it can step
        }
        // at the last (100-th) step, it returns false to stop the simulation.
        BOOST_TEST(!sim->step());
        sim->finalize();
    }
}

BOOST_AUTO_TEST_CASE(read_langevin_molecular_dynamics_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_simulator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    auto root = mjolnir::test::make_empty_input();

    {
        using namespace toml::literals;
        const auto v = u8R"(
            type            = "MolecularDynamics"
            integrator.type = "UnderdampedLangevin"
            integrator.seed = 12345
            integrator.parameters = []
            precision       = "double"
            boundary_type   = "Unlimited"
            delta_t         = 0.1
            total_step      = 100
            save_step       = 10
        )"_toml;

        root["simulator"] = v;
        const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto mdsim = dynamic_cast<mjolnir::MolecularDynamicsSimulator<
            traits_type, mjolnir::UnderdampedLangevinIntegrator<traits_type>>*>(sim.get());
        BOOST_TEST(static_cast<bool>(mdsim));

        sim->initialize();
        for(std::size_t i=0; i<99; ++i)
        {
            BOOST_TEST(mdsim->time() == i * 0.1, boost::test_tools::tolerance(tol));
            BOOST_TEST(sim->step()); // check it can step
        }
        // at the last (100-th) step, it returns false to stop the simulation.
        BOOST_TEST(!sim->step());
        sim->finalize();
    }
}

BOOST_AUTO_TEST_CASE(read_BAOAB_langevin_molecular_dynamics_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_simulator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    auto root = mjolnir::test::make_empty_input();

    {
        using namespace toml::literals;
        const auto v = u8R"(
            type            = "MolecularDynamics"
            integrator.type = "BAOABLangevin"
            integrator.seed = 12345
            integrator.parameters = []
            precision       = "double"
            boundary_type   = "Unlimited"
            delta_t         = 0.1
            total_step      = 100
            save_step       = 10
        )"_toml;

        root["simulator"] = v;
        const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto mdsim = dynamic_cast<mjolnir::MolecularDynamicsSimulator<
            traits_type, mjolnir::BAOABLangevinIntegrator<traits_type>>*>(sim.get());
        BOOST_TEST(static_cast<bool>(mdsim));

        sim->initialize();
        for(std::size_t i=0; i<99; ++i)
        {
            BOOST_TEST(mdsim->time() == i * 0.1, boost::test_tools::tolerance(tol));
            BOOST_TEST(sim->step()); // check it can step
        }
        // at the last (100-th) step, it returns false to stop the simulation.
        BOOST_TEST(!sim->step());
        sim->finalize();
    }
}

BOOST_AUTO_TEST_CASE(read_steepest_descent_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_simulator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    auto root = mjolnir::test::make_empty_input();
    {
        using namespace toml::literals;
        const auto v = u8R"(
            type            = "SteepestDescent"
            precision       = "double"
            boundary_type   = "Unlimited"
            delta           = 0.1
            step_limit      = 100
            save_step       = 10
            threshold       = 0.0
        )"_toml;

        root["simulator"] = v;
        const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto sdsim = dynamic_cast<
            mjolnir::SteepestDescentSimulator<traits_type>*>(sim.get());
        BOOST_TEST(static_cast<bool>(sdsim));

        sim->initialize();
        for(std::size_t i=0; i<99; ++i)
        {
            BOOST_TEST(sim->step()); // check it can step
        }
        // at the last (100-th) step, it returns false to stop the simulation.
        BOOST_TEST(!sim->step());
        sim->finalize();
    }
}

BOOST_AUTO_TEST_CASE(read_simulated_annealing_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_simulator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    auto root = mjolnir::test::make_empty_input();
    {
        using namespace toml::literals;
        const auto v = u8R"(
            type            = "SimulatedAnnealing"
            integrator.type = "UnderdampedLangevin"
            integrator.seed = 12345
            integrator.parameters = []
            precision       = "double"
            boundary_type   = "Unlimited"
            total_step      = 100
            save_step       = 10
            each_step       = 1
            delta_t         = 0.1
            schedule.type  = "linear"
            schedule.begin = 300.0
            schedule.end   =  10.0

        )"_toml;

        root["simulator"] = v;
        const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto sasim = dynamic_cast<mjolnir::SimulatedAnnealingSimulator<
            traits_type, mjolnir::UnderdampedLangevinIntegrator<traits_type>,
            mjolnir::LinearScheduler>*>(sim.get());
        BOOST_TEST(static_cast<bool>(sasim));

        sim->initialize();
        for(std::size_t i=0; i<99; ++i)
        {
            BOOST_TEST(sasim->system().attribute("temperature") ==
                       300.0 * ((100-i) / 100.0) + 10.0 * (i / 100.0),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(sim->step()); // check it can step
        }
        // at the last (100-th) step, it returns false to stop the simulation.
        BOOST_TEST(!sim->step());
        sim->finalize();
    }
}

BOOST_AUTO_TEST_CASE(read_BAOAB_Langevin_simulated_annealing_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_simulator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    auto root = mjolnir::test::make_empty_input();
    {
        using namespace toml::literals;
        const auto v = u8R"(
            type            = "SimulatedAnnealing"
            integrator.type = "BAOABLangevin"
            integrator.seed = 12345
            integrator.parameters = []
            precision       = "double"
            boundary_type   = "Unlimited"
            total_step      = 100
            save_step       = 10
            each_step       = 1
            delta_t         = 0.1
            schedule.type  = "linear"
            schedule.begin = 300.0
            schedule.end   =  10.0

        )"_toml;

        root["simulator"] = v;
        const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto sasim = dynamic_cast<mjolnir::SimulatedAnnealingSimulator<
            traits_type, mjolnir::BAOABLangevinIntegrator<traits_type>,
            mjolnir::LinearScheduler>*>(sim.get());
        BOOST_TEST(static_cast<bool>(sasim));

        sim->initialize();
        for(std::size_t i=0; i<99; ++i)
        {
            BOOST_TEST(sasim->system().attribute("temperature") ==
                       300.0 * ((100-i) / 100.0) + 10.0 * (i / 100.0),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(sim->step()); // check it can step
        }
        // at the last (100-th) step, it returns false to stop the simulation.
        BOOST_TEST(!sim->step());
        sim->finalize();
    }
}
