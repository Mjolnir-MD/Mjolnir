#define BOOST_TEST_MODULE "test_read_input_path"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_files_table.hpp>

BOOST_AUTO_TEST_CASE(read_input_path)
{
    mjolnir::LoggerManager::set_default_logger("test_read_input_path.log");
    {
        const toml::table t = toml::table{
            {"files", toml::table{
                {"input", toml::table{
                    {"path", "./input"}
                }}
            }}
        };
        const auto input_path = mjolnir::read_input_path(t);
        BOOST_TEST("./input/");
    }
    {
        const toml::table t = toml::table{
            {"files", toml::table{
                {"input", toml::table{
                    {"path", "./input/"}
                }}
            }}
        };
        const auto input_path = mjolnir::read_input_path(t);
        BOOST_TEST("./input/");
    }
    {
        const toml::table t = toml::table{
            {"files", toml::table{
                {"input", toml::table{}}
            }}
        };
        const auto input_path = mjolnir::read_input_path(t);
        BOOST_TEST("./");
    }
    {
        const toml::table t = toml::table{
            {"files", toml::table{}}
        };
        const auto input_path = mjolnir::read_input_path(t);
        BOOST_TEST("./");
    }
}
