#define BOOST_TEST_MODULE "test_read_steepest_descent_simulator"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <test/util/make_empty_input.hpp>
#include <mjolnir/input/read_simulator.hpp>

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
