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
        const toml::table v = toml::table{{"files", toml::table{{"output",
                toml::table{{"path", "./"}, {"prefix", "test"}, {"format", "xyz"}}
            }}
        }};

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);
        BOOST_TEST(obs.observers().at(0)->prefix() == "./test");
        BOOST_TEST(obs.observers().at(1)->prefix() == "./test");

        const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(
                             obs.observers().at(0).get());
        BOOST_TEST(static_cast<bool>(ene));

        const auto xyz = dynamic_cast<mjolnir::XYZObserver<traits_type>*>(
                             obs.observers().at(1).get());
        BOOST_TEST(static_cast<bool>(xyz));
    }
    {
        const toml::table v = toml::table{{"files", toml::table{{"output",
                toml::table{{"prefix", "test"}, {"format", "xyz"}}
            }}
        }};

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);
        BOOST_TEST(obs.observers().at(0)->prefix() == "./test");
        BOOST_TEST(obs.observers().at(1)->prefix() == "./test");

        const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(
                             obs.observers().at(0).get());
        BOOST_TEST(static_cast<bool>(ene));

        const auto xyz = dynamic_cast<mjolnir::XYZObserver<traits_type>*>(
                             obs.observers().at(1).get());
        BOOST_TEST(static_cast<bool>(xyz));
    }

    // XXX
    // The following block tests read_observer successfully adds a path before
    // the prefix. But it requires the directory `test` under WORKING_DIRECTORY.
    // It strongly depends on the directory structure...
    // It is not good. We need to find a way to avoid this dependency.
    {
        const toml::table v = toml::table{{"files", toml::table{{"output",
                toml::table{{"path", "./test"}, {"prefix", "test"}, {"format", "xyz"}}
            }}
        }};

        const auto obs = mjolnir::read_observer<traits_type>(v);

        BOOST_TEST(obs.observers().size() == 2u);
        BOOST_TEST(obs.observers().at(0)->prefix() == "./test/test");
        BOOST_TEST(obs.observers().at(1)->prefix() == "./test/test");

        const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(
                             obs.observers().at(0).get());
        BOOST_TEST(static_cast<bool>(ene));

        const auto xyz = dynamic_cast<mjolnir::XYZObserver<traits_type>*>(
                             obs.observers().at(1).get());
        BOOST_TEST(static_cast<bool>(xyz));
    }
}

// check that read_observer returns XYZObserver or not
BOOST_AUTO_TEST_CASE(read_xyz_observer)
{
    mjolnir::LoggerManager::set_default_logger("test_read_observer.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        const toml::table files{
            {"output", toml::table{{"path", "./"}, {"prefix", "test"}, {"format", "xyz"}}}
        };

        const toml::table v{{"files", files}};

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);

        const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(
                             obs.observers().at(0).get());
        BOOST_TEST(static_cast<bool>(ene));

        const auto xyz = dynamic_cast<mjolnir::XYZObserver<traits_type>*>(
                             obs.observers().at(1).get());
        BOOST_TEST(static_cast<bool>(xyz));
    }
}

// check that read_observer returns DCDObserver or not
BOOST_AUTO_TEST_CASE(read_dcd_observer)
{
    mjolnir::LoggerManager::set_default_logger("test_read_observer.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        const toml::table files{
            {"output", toml::table{{"path", "./"}, {"prefix", "test"}, {"format", "dcd"}}}
        };

        const toml::table v{{"files", files}};

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);

        const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(
                             obs.observers().at(0).get());
        BOOST_TEST(static_cast<bool>(ene));

        const auto dcd = dynamic_cast<mjolnir::DCDObserver<traits_type>*>(
                             obs.observers().at(1).get());
        BOOST_TEST(static_cast<bool>(dcd));
    }
}

// check that read_observer returns TRRObserver or not
BOOST_AUTO_TEST_CASE(read_trr_observer)
{
    mjolnir::LoggerManager::set_default_logger("test_read_observer.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        const toml::table files{
            {"output", toml::table{{"path", "./"}, {"prefix", "test"}, {"format", "trr"}}}
        };

        const toml::table v{{"files", files}};

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.observers().size() == 2u);

        const auto ene = dynamic_cast<mjolnir::EnergyObserver<traits_type>*>(
                             obs.observers().at(0).get());
        BOOST_TEST(static_cast<bool>(ene));

        const auto trr = dynamic_cast<mjolnir::TRRObserver<traits_type>*>(
                             obs.observers().at(1).get());
        BOOST_TEST(static_cast<bool>(trr));
    }
}
