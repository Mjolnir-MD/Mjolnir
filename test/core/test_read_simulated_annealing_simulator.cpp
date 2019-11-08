#define BOOST_TEST_MODULE "test_read_simulated_annealing_simulator"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <test/util/make_empty_input.hpp>
#include <mjolnir/input/read_simulator.hpp>

BOOST_AUTO_TEST_CASE(read_simulated_annealing_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_simulated_annealing_simulator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using integrator_type = mjolnir::UnderdampedLangevinIntegrator<traits_type>;

    constexpr real_type tol = 1e-8;
    auto root = mjolnir::test::make_empty_input();

    using namespace toml::literals;
    const auto v = u8R"(
        type            = "SimulatedAnnealing"
        integrator.type = "UnderdampedLangevin"
        integrator.parameters = []
        precision       = "double"
        boundary_type   = "Unlimited"
        seed            = 12345
        total_step      = 100
        save_step       = 10
        each_step       = 1
        delta_t         = 0.1
        schedule.type  = "linear"
        schedule.begin = 300.0
        schedule.end   =  10.0

    )"_toml;
    root.as_table()["simulator"] = v;

    {
        const auto sim = mjolnir::read_simulator<traits_type, integrator_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto sasim = dynamic_cast<mjolnir::SimulatedAnnealingSimulator<
            traits_type, integrator_type, mjolnir::LinearScheduler>*>(sim.get());
        BOOST_TEST(static_cast<bool>(sasim));

        BOOST_TEST(sasim->rng().seed() == 12345u);

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
    {
        const auto sim = mjolnir::read_integrator_type<traits_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto sasim = dynamic_cast<mjolnir::SimulatedAnnealingSimulator<
            traits_type, integrator_type, mjolnir::LinearScheduler>*>(sim.get());
        BOOST_TEST(static_cast<bool>(sasim));

        BOOST_TEST(sasim->rng().seed() == 12345u);

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
    mjolnir::LoggerManager::set_default_logger("test_read_simulated_annealing_simulator.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using integrator_type = mjolnir::BAOABLangevinIntegrator<traits_type>;

    constexpr real_type tol = 1e-8;
    auto root = mjolnir::test::make_empty_input();

    using namespace toml::literals;
    const auto v = u8R"(
        type            = "SimulatedAnnealing"
        integrator.type = "BAOABLangevin"
        integrator.parameters = []
        precision       = "double"
        boundary_type   = "Unlimited"
        seed            = 12345
        total_step      = 100
        save_step       = 10
        each_step       = 1
        delta_t         = 0.1
        schedule.type  = "linear"
        schedule.begin = 300.0
        schedule.end   =  10.0

    )"_toml;

    root.as_table()["simulator"] = v;
    {
        const auto sim = mjolnir::read_simulator<traits_type, integrator_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto sasim = dynamic_cast<mjolnir::SimulatedAnnealingSimulator<
            traits_type, integrator_type, mjolnir::LinearScheduler>*>(sim.get());
        BOOST_TEST(static_cast<bool>(sasim));

        BOOST_TEST(sasim->rng().seed() == 12345u);

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
    {
        const auto sim = mjolnir::read_integrator_type<traits_type>(root, v);
        BOOST_TEST(static_cast<bool>(sim));

        const auto sasim = dynamic_cast<mjolnir::SimulatedAnnealingSimulator<
            traits_type, integrator_type, mjolnir::LinearScheduler>*>(sim.get());
        BOOST_TEST(static_cast<bool>(sasim));

        BOOST_TEST(sasim->rng().seed() == 12345u);

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
