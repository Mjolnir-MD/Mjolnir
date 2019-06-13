#define BOOST_TEST_MODULE "test_read_path"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_path.hpp>

BOOST_AUTO_TEST_CASE(read_input_path)
{
    mjolnir::LoggerManager::set_default_logger("test_read_input_path.log");
    {
        using namespace toml::literals;
        const toml::table t = toml::get<toml::table>(u8R"(
            files.input.path = "./input"
        )"_toml);

        const auto input_path = mjolnir::read_input_path(t);
        BOOST_TEST(input_path == "./input/");
    }
    {
       using namespace toml::literals;
        const toml::table t = toml::get<toml::table>(u8R"(
            files.input.path = "./input/"
        )"_toml);

        const auto input_path = mjolnir::read_input_path(t);
        BOOST_TEST(input_path == "./input/");
    }
    {
       using namespace toml::literals;
        const toml::table t = toml::get<toml::table>(u8R"(
            files.input = {}
        )"_toml);

        const auto input_path = mjolnir::read_input_path(t);
        BOOST_TEST(input_path == "./");
    }
    {
       using namespace toml::literals;
        const toml::table t = toml::get<toml::table>(u8R"(
            files = {}
        )"_toml);

        const auto input_path = mjolnir::read_input_path(t);
        BOOST_TEST(input_path == "./");
    }
}

BOOST_AUTO_TEST_CASE(read_output_path)
{
    mjolnir::LoggerManager::set_default_logger("test_read_output_path.log");
    {
        using namespace toml::literals;
        const toml::table t = toml::get<toml::table>(u8R"(
            files.output.path = "./output"
        )"_toml);

        const auto output_path = mjolnir::read_output_path(t);
        BOOST_TEST(output_path == "./output/");
    }
    {
       using namespace toml::literals;
        const toml::table t = toml::get<toml::table>(u8R"(
            files.output.path = "./output/"
        )"_toml);

        const auto output_path = mjolnir::read_output_path(t);
        BOOST_TEST(output_path == "./output/");
    }
    {
       using namespace toml::literals;
        const toml::table t = toml::get<toml::table>(u8R"(
            files.output = {}
        )"_toml);

        const auto output_path = mjolnir::read_output_path(t);
        BOOST_TEST(output_path == "./");
    }
    {
       using namespace toml::literals;
        const toml::table t = toml::get<toml::table>(u8R"(
            files = {}
        )"_toml);

        const auto output_path = mjolnir::read_output_path(t);
        BOOST_TEST(output_path == "./");
    }
}

