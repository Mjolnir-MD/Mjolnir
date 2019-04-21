#define BOOST_TEST_MODULE "test_bond_angle_interaction"
#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_utility.hpp>

BOOST_AUTO_TEST_CASE(test_read_array_from_table)
{
    mjolnir::LoggerManager::set_default_logger("test_read_utility.log");
    const auto reader = [](const toml::value& v) -> std::pair<double, int> {
        const auto a = toml::find<double>(v, "a");
        const auto b = toml::find<int   >(v, "b");
        return std::make_pair(a, b);
    };

    using namespace toml::literals;
    {
        const toml::value v = u8R"(
            size    = 10
            default = {a = 1.0, b = 10}
            )"_toml;

        const auto result =
            mjolnir::read_array<std::pair<double, int>>(v, reader);

        BOOST_TEST(result.size()  == 10u);
        for(const auto& item : result)
        {
            BOOST_TEST(item.first  == 1.0);
            BOOST_TEST(item.second ==  10);
        }
    }
    {
        const toml::value v = u8R"(
            size    = 10
            default = {a = 1.0, b = 10}
            values  = [
                # overwrite special values
                {index = 0, a = 3.14, b = 42  },
                {index = 5, b = 54,   a = 2.71},
            ]
            )"_toml;

        const auto result =
            mjolnir::read_array<std::pair<double, int>>(v, reader);

        BOOST_TEST(result.size()  == 10u);
        for(std::size_t i=0; i<10; ++i)
        {
            const auto& item = result.at(i);
            if(i == 0)
            {
                BOOST_TEST(item.first  == 3.14);
                BOOST_TEST(item.second ==   42);
            }
            else if(i == 5)
            {
                BOOST_TEST(item.first  == 2.71);
                BOOST_TEST(item.second ==   54);
            }
            else
            {
                BOOST_TEST(item.first  == 1.0);
                BOOST_TEST(item.second ==  10);
            }
        }
    }
    {
        const toml::value v = u8R"(
            size = 10
            # no default value
            )"_toml;

        bool thrown = false;
        try
        {
            const auto result =
                mjolnir::read_array<std::pair<double, int>>(v, reader);
        }
        catch(const std::runtime_error& re)
        {
            std::cout << "what: " << re.what() << std::endl;
            thrown = true;
        }
        BOOST_TEST(thrown);
    }
    {
        const toml::value v = u8R"(
            # no size
            default = {a = 1.0, b = 10}
            )"_toml;

        bool thrown = false;
        try
        {
            const auto result =
                mjolnir::read_array<std::pair<double, int>>(v, reader);
        }
        catch(const std::runtime_error& re)
        {
            std::cout << "what: " << re.what() << std::endl;
            thrown = true;
        }
        BOOST_TEST(thrown);
    }
    {
        const toml::value v = u8R"(
            size    = 10
            default = {a = 1.0, b = 10}
            values  = [
                {index = 100, a = 10.0, b = 42} # index exceeds the size
            ]
            )"_toml;

        bool thrown = false;
        try
        {
            const auto result =
                mjolnir::read_array<std::pair<double, int>>(v, reader);
        }
        catch(const std::runtime_error& re)
        {
            std::cout << "what: " << re.what() << std::endl;
            thrown = true;
        }
        BOOST_TEST(thrown);
    }
}

BOOST_AUTO_TEST_CASE(test_read_array_from_array)
{
    mjolnir::LoggerManager::set_default_logger("test_read_utility.log");
    const auto reader = [](const toml::value& v) -> std::pair<double, int> {
        const auto a = toml::find<double>(v, "a");
        const auto b = toml::find<int   >(v, "b");
        return std::make_pair(a, b);
    };

    using namespace toml::literals;
    {
        const toml::value v = u8R"([
                {index = 0, a = 1.00, b = 10},
                {index = 1, a = 3.14, b = 42},
                {index = 2, a = 2.71, b = 54},
            ])"_toml;

        const auto result =
            mjolnir::read_array<std::pair<double, int>>(v, reader);

        BOOST_TEST(result.size()  == 3u);
        BOOST_TEST(result.at(0).first     == 1.00);
        BOOST_TEST(result.at(0).second    ==   10);
        BOOST_TEST(result.at(1).first     == 3.14);
        BOOST_TEST(result.at(1).second    ==   42);
        BOOST_TEST(result.at(2).first     == 2.71);
        BOOST_TEST(result.at(2).second    ==   54);
    }

    {
        const toml::value v = u8R"([
                {index = 100, a = 1.00, b = 10}, # range exceeds
            ])"_toml;

        bool thrown = false;
        try
        {
            const auto result =
                mjolnir::read_array<std::pair<double, int>>(v, reader);
        }
        catch(const std::runtime_error& re)
        {
            std::cout << "what: " << re.what() << std::endl;
            thrown = true;
        }
        BOOST_TEST(thrown);
    }

    {
        const toml::value v = u8R"([
                {index = 0, a = 1.00, b = 10}, # range [0, 1] not exhaustive
                {index = 0, a = 1.00, b = 10}, # because it is also 0-th value
            ])"_toml;

        bool thrown = false;
        try
        {
            const auto result =
                mjolnir::read_array<std::pair<double, int>>(v, reader);
        }
        catch(const std::runtime_error& re)
        {
            std::cout << "what: " << re.what() << std::endl;
            thrown = true;
        }
        BOOST_TEST(thrown);
    }
}
