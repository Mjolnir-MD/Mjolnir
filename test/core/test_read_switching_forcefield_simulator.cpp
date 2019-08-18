#define BOOST_TEST_MODULE "test_read_switching_forcefield_simulator"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <test/util/make_empty_input.hpp>
#include <mjolnir/input/read_simulator.hpp>

BOOST_AUTO_TEST_CASE(read_newtonian_switching_forcefield_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_switching_forcefield_simulator.log");

    using real_type   = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    auto root = mjolnir::test::make_empty_input();

    using namespace toml::literals;

    const auto forcefields = toml::find(u8R"(
        [[forcefields]]
        name = "open"
        [[forcefields]]
        name = "close"
        [[forcefields.local]]
        interaction = "BondLength"
        potential = "Harmonic"
        topology  = "bond"
        parameters = [
            {indices = [0, 1], v0 = 1.0, k = 1.0}
        ]
    )"_toml, "forcefields");
    root.as_table()["forcefields"] = forcefields;

    const auto simulator = u8R"(
        type            = "SwitchingForceField"
        integrator.type = "VelocityVerlet"
        precision       = "double"
        boundary_type   = "Unlimited"
        delta_t         = 0.1
        total_step      = 30
        save_step       =  1
        schedule        = [
            {until = 10, forcefield = "open"},
            {until = 20, forcefield = "close"},
            {until = 30, forcefield = "open"},
        ]
    )"_toml;
    root.as_table()["simulator"] = simulator;

    const auto system = u8R"(
        attributes.temperature = 300.0
        boundary_shape = {} # Unlimited
        particles = [
            {mass=1.0, position=[0.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]},
            {mass=1.0, position=[2.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]},
        ]
    )"_toml;
    root.as_table()["systems"] = toml::array{system};

    const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, simulator);
    BOOST_TEST(static_cast<bool>(sim));

    const auto mdsim = dynamic_cast<mjolnir::SwitchingForceFieldSimulator<
        traits_type, mjolnir::VelocityVerletIntegrator<traits_type>>*>(sim.get());
    BOOST_TEST(static_cast<bool>(mdsim));

    BOOST_TEST(mdsim->forcefield_index().at("open" ) == 0u);
    BOOST_TEST(mdsim->forcefield_index().at("close") == 1u);

    for(const auto& kv : mdsim->schedule())
    {
        std::cout << '{' << kv.first << ':' << kv.second << '}' << std::endl;
    }

    BOOST_TEST(mdsim->schedule().size() == 3u);
    BOOST_TEST(mdsim->schedule().at(0).first  == 10u);
    BOOST_TEST(mdsim->schedule().at(0).second == "open");
    BOOST_TEST(mdsim->schedule().at(1).first  == 20u);
    BOOST_TEST(mdsim->schedule().at(1).second == "close");
    BOOST_TEST(mdsim->schedule().at(2).first  == 30u);
    BOOST_TEST(mdsim->schedule().at(2).second == "open");

    BOOST_TEST(mdsim->system().size() == 2u);

    mdsim->initialize();
    for(std::size_t i=0; i<10; ++i)
    {
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) == 0.0);
        mdsim->step();
    }
    for(std::size_t i=0; i<10; ++i)
    {
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) != 0.0);
        mdsim->step();
    }
    for(std::size_t i=0; i<10; ++i)
    {
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) == 0.0);
        mdsim->step();
    }
    mdsim->finalize();
}

BOOST_AUTO_TEST_CASE(read_underdamped_langevin_switching_forcefield_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_switching_forcefield_simulator.log");

    using real_type   = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    auto root = mjolnir::test::make_empty_input();

    using namespace toml::literals;

    const auto forcefields = toml::find(u8R"(
        [[forcefields]]
        name = "open"
        [[forcefields]]
        name = "close"
        [[forcefields.local]]
        interaction = "BondLength"
        potential = "Harmonic"
        topology  = "bond"
        parameters = [
            {indices = [0, 1], v0 = 1.0, k = 1.0}
        ]
    )"_toml, "forcefields");
    root.as_table()["forcefields"] = forcefields;

    const auto simulator = u8R"(
        type            = "SwitchingForceField"
        integrator.type = "UnderdampedLangevin"
        integrator.seed = 1
        integrator.parameters = [{index = 0, gamma = 1.0}, {index = 1, gamma = 1.0}]
        precision       = "double"
        boundary_type   = "Unlimited"
        delta_t         = 0.1
        total_step      = 30
        save_step       =  1
        schedule        = [
            {until = 10, forcefield = "open"},
            {until = 20, forcefield = "close"},
            {until = 30, forcefield = "open"},
        ]
    )"_toml;
    root.as_table()["simulator"] = simulator;

    const auto system = u8R"(
        # set temperature zero to make random force zero
        attributes.temperature = 0.0
        boundary_shape = {} # Unlimited
        particles = [
            {mass=1.0, position=[0.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]},
            {mass=1.0, position=[2.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]},
        ]
    )"_toml;
    root.as_table()["systems"] = toml::array{system};

    const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, simulator);
    BOOST_TEST(static_cast<bool>(sim));

    const auto mdsim = dynamic_cast<mjolnir::SwitchingForceFieldSimulator<
        traits_type, mjolnir::UnderdampedLangevinIntegrator<traits_type>>*>(sim.get());
    BOOST_TEST(static_cast<bool>(mdsim));

    BOOST_TEST(mdsim->forcefield_index().at("open" ) == 0u);
    BOOST_TEST(mdsim->forcefield_index().at("close") == 1u);

    for(const auto& kv : mdsim->schedule())
    {
        std::cout << '{' << kv.first << ':' << kv.second << '}' << std::endl;
    }

    BOOST_TEST(mdsim->schedule().size() == 3u);
    BOOST_TEST(mdsim->schedule().at(0).first  == 10u);
    BOOST_TEST(mdsim->schedule().at(0).second == "open");
    BOOST_TEST(mdsim->schedule().at(1).first  == 20u);
    BOOST_TEST(mdsim->schedule().at(1).second == "close");
    BOOST_TEST(mdsim->schedule().at(2).first  == 30u);
    BOOST_TEST(mdsim->schedule().at(2).second == "open");

    BOOST_TEST(mdsim->system().size() == 2u);

    mdsim->initialize();
    for(std::size_t i=0; i<10; ++i)
    {
        // empty forcefield. no forces.
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) == 0.0);
        mdsim->step();
    }
    for(std::size_t i=0; i<10; ++i)
    {
        // harmonic potential is applied. non-zero force.
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) != 0.0);
        mdsim->step();
    }
    for(std::size_t i=0; i<10; ++i)
    {
        // empty forcefield. no forces.
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) == 0.0);
        mdsim->step();
    }
    mdsim->finalize();
}

BOOST_AUTO_TEST_CASE(read_BAOAB_langevin_switching_forcefield_simulator)
{
    mjolnir::LoggerManager::set_default_logger("test_read_switching_forcefield_simulator.log");

    using real_type   = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    auto root = mjolnir::test::make_empty_input();

    using namespace toml::literals;

    const auto forcefields = toml::find(u8R"(
        [[forcefields]]
        name = "open"
        [[forcefields]]
        name = "close"
        [[forcefields.local]]
        interaction = "BondLength"
        potential = "Harmonic"
        topology  = "bond"
        parameters = [
            {indices = [0, 1], v0 = 1.0, k = 1.0}
        ]
    )"_toml, "forcefields");
    root.as_table()["forcefields"] = forcefields;

    const auto simulator = u8R"(
        type            = "SwitchingForceField"
        integrator.type = "BAOABLangevin"
        integrator.seed = 1
        integrator.parameters = [{index = 0, gamma = 1.0}, {index = 1, gamma = 1.0}]
        precision       = "double"
        boundary_type   = "Unlimited"
        delta_t         = 0.1
        total_step      = 30
        save_step       =  1
        schedule        = [
            {until = 10, forcefield = "open"},
            {until = 20, forcefield = "close"},
            {until = 30, forcefield = "open"},
        ]
    )"_toml;
    root.as_table()["simulator"] = simulator;

    const auto system = u8R"(
        # set temperature zero to make random force zero
        attributes.temperature = 0.0
        boundary_shape = {} # Unlimited
        particles = [
            {mass=1.0, position=[0.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]},
            {mass=1.0, position=[2.0, 0.0, 0.0], velocity=[0.0, 0.0, 0.0]},
        ]
    )"_toml;
    root.as_table()["systems"] = toml::array{system};

    const auto sim = mjolnir::read_simulator_from_table<traits_type>(root, simulator);
    BOOST_TEST(static_cast<bool>(sim));

    const auto mdsim = dynamic_cast<mjolnir::SwitchingForceFieldSimulator<
        traits_type, mjolnir::BAOABLangevinIntegrator<traits_type>>*>(sim.get());
    BOOST_TEST(static_cast<bool>(mdsim));

    BOOST_TEST(mdsim->forcefield_index().at("open" ) == 0u);
    BOOST_TEST(mdsim->forcefield_index().at("close") == 1u);

    for(const auto& kv : mdsim->schedule())
    {
        std::cout << '{' << kv.first << ':' << kv.second << '}' << std::endl;
    }

    BOOST_TEST(mdsim->schedule().size() == 3u);
    BOOST_TEST(mdsim->schedule().at(0).first  == 10u);
    BOOST_TEST(mdsim->schedule().at(0).second == "open");
    BOOST_TEST(mdsim->schedule().at(1).first  == 20u);
    BOOST_TEST(mdsim->schedule().at(1).second == "close");
    BOOST_TEST(mdsim->schedule().at(2).first  == 30u);
    BOOST_TEST(mdsim->schedule().at(2).second == "open");

    BOOST_TEST(mdsim->system().size() == 2u);

    mdsim->initialize();
    for(std::size_t i=0; i<10; ++i)
    {
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) == 0.0);
        mdsim->step();
    }
    for(std::size_t i=0; i<10; ++i)
    {
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) != 0.0);
        mdsim->step();
    }
    for(std::size_t i=0; i<10; ++i)
    {
        BOOST_TEST(mjolnir::math::X(mdsim->system().force(0)) == 0.0);
        mdsim->step();
    }
    mdsim->finalize();
}
