#define BOOST_TEST_MODULE "test_file_inclusion"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/utility.hpp>

BOOST_AUTO_TEST_CASE(test_file_inclusion_normal)
{
    mjolnir::LoggerManager::set_default_logger("test_file_inclusion.log");
    toml::value includes = toml::table{
        {"include", "included.toml"}
    };
    toml::value included = toml::table{
        {"foo", 3.14}, {"bar", 42}
    };
    {
        std::ofstream ofs("included.toml");
        ofs << included;
    }

    mjolnir::expand_include(includes);

    BOOST_TEST(toml::find<double>(includes, "foo") == std::stod("3.14"));
    BOOST_TEST(toml::find<int   >(includes, "bar") == 42);
}
BOOST_AUTO_TEST_CASE(test_file_inclusion_array_of_tables)
{
    mjolnir::LoggerManager::set_default_logger("test_file_inclusion.log");
    toml::value includes = toml::array{
        toml::table{{"include", "included.toml"}}
    };
    toml::value included = toml::table{
        {"foo", 3.14}, {"bar", 42}
    };
    {
        std::ofstream ofs("included.toml");
        ofs << included;
    }

    mjolnir::expand_include(includes);

    BOOST_TEST(toml::find<double>(includes.at(0), "foo") == std::stod("3.14"));
    BOOST_TEST(toml::find<int   >(includes.at(0), "bar") == 42);
}
