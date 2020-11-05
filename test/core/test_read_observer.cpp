#define BOOST_TEST_MODULE "test_read_observer"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_observer.hpp>

BOOST_AUTO_TEST_CASE(read_observer)
{
    mjolnir::LoggerManager::set_default_logger("test_read_observer.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            [files]
            output.path   = "./"
            output.prefix = "test"
            output.format = "xyz"
        )"_toml);

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);
        BOOST_TEST(obs.observers().at(0)->prefix() == "./test");
        BOOST_TEST(obs.observers().at(1)->prefix() == "./test");

        bool has_xyz = false;
        bool has_ene = false;
        for(const auto& ptr : obs.observers())
        {
            const auto xyz = dynamic_cast<mjolnir::XYZObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(xyz))
            {
                has_xyz = true;
            }
            const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(ene))
            {
                has_ene = true;
            }
        }

        BOOST_TEST(has_xyz);
        BOOST_TEST(has_ene);
    }
    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            [files]
            output.prefix = "test"
            output.format = "xyz"
        )"_toml);

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);
        BOOST_TEST(obs.observers().at(0)->prefix() == "./test");
        BOOST_TEST(obs.observers().at(1)->prefix() == "./test");

        bool has_xyz = false;
        bool has_ene = false;
        for(const auto& ptr : obs.observers())
        {
            const auto xyz = dynamic_cast<mjolnir::XYZObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(xyz))
            {
                has_xyz = true;
            }
            const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(ene))
            {
                has_ene = true;
            }
        }
        BOOST_TEST(has_xyz);
        BOOST_TEST(has_ene);
    }

    // XXX
    // The following block tests read_observer successfully adds a path before
    // the prefix. But it requires the directory `test` under WORKING_DIRECTORY.
    // It strongly depends on the directory structure...
    // It is not good. We need to find a way to avoid this dependency.
    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            [files]
            output.path   = "./test"
            output.prefix = "test"
            output.format = "xyz"
        )"_toml);

        const auto obs = mjolnir::read_observer<traits_type>(v);

        BOOST_TEST(obs.observers().size() == 2u);
        BOOST_TEST(obs.observers().at(0)->prefix() == "./test/test");
        BOOST_TEST(obs.observers().at(1)->prefix() == "./test/test");

        bool has_xyz = false;
        bool has_ene = false;
        for(const auto& ptr : obs.observers())
        {
            const auto xyz = dynamic_cast<mjolnir::XYZObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(xyz))
            {
                has_xyz = true;
            }
            const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(ene))
            {
                has_ene = true;
            }
        }
        BOOST_TEST(has_xyz);
        BOOST_TEST(has_ene);
    }
}

// check that read_observer returns XYZObserver or not
BOOST_AUTO_TEST_CASE(read_xyz_observer)
{
    mjolnir::LoggerManager::set_default_logger("test_read_observer.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            [files]
            output.path   = "./"
            output.prefix = "test"
            output.format = "xyz"
        )"_toml);

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);

        bool has_xyz = false;
        bool has_ene = false;
        for(const auto& ptr : obs.observers())
        {
            const auto xyz = dynamic_cast<mjolnir::XYZObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(xyz))
            {
                has_xyz = true;
            }
            const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(ene))
            {
                has_ene = true;
            }
        }
        BOOST_TEST(has_xyz);
        BOOST_TEST(has_ene);
    }
}

// check that read_observer returns DCDObserver or not
BOOST_AUTO_TEST_CASE(read_dcd_observer)
{
    mjolnir::LoggerManager::set_default_logger("test_read_observer.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            [files]
            output.path   = "./"
            output.prefix = "test"
            output.format = "dcd"
        )"_toml);

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);

        bool has_dcd = false;
        bool has_ene = false;
        for(const auto& ptr : obs.observers())
        {
            const auto dcd = dynamic_cast<mjolnir::DCDObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(dcd))
            {
                has_dcd = true;
            }
            const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(ene))
            {
                has_ene = true;
            }
        }
        BOOST_TEST(has_dcd);
        BOOST_TEST(has_ene);
    }
}

// check that read_observer returns TRRObserver or not
BOOST_AUTO_TEST_CASE(read_trr_observer)
{
    mjolnir::LoggerManager::set_default_logger("test_read_observer.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            [files]
            output.path   = "./"
            output.prefix = "test"
            output.format = "trr"
        )"_toml);

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);

        bool has_trr = false;
        bool has_ene = false;
        for(const auto& ptr : obs.observers())
        {
            const auto trr = dynamic_cast<mjolnir::TRRObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(trr))
            {
                has_trr = true;
            }
            const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(ptr.get());
            if(static_cast<bool>(ene))
            {
                has_ene = true;
            }
        }
        BOOST_TEST(has_trr);
        BOOST_TEST(has_ene);
    }
}
