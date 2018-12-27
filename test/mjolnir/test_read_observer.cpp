#define BOOST_TEST_MODULE "test_read_observer"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_observer.hpp>

BOOST_AUTO_TEST_CASE(read_observer)
{
    mjolnir::LoggerManager::set_default_logger("test_read_observer.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{{"files", toml::table{{"output",
                toml::table{{"path", "./"}, {"prefix", "test"}}
            }}
        }};

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.prefix() == "./test");
    }
    {
        const toml::table v = toml::table{{"files", toml::table{{"output",
                toml::table{{"prefix", "test"}}
            }}
        }};

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.prefix() == "./test");
    }

    // XXX
    // The following block tests read_observer successfully adds a path before
    // the prefix. But it requires the directory `test` under WORKING_DIRECTORY.
    // It strongly depends on the directory structure...
    // It is not good. We need to find a way to avoid this dependency.
    {
        const toml::table v = toml::table{{"files", toml::table{{"output",
                toml::table{{"path", "./test"}, {"prefix", "test"}}
            }}
        }};

        const auto obs = mjolnir::read_observer<traits_type>(v);
        BOOST_TEST(obs.prefix() == "./test/test");
    }
}