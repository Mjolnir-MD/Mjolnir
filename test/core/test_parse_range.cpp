#define BOOST_TEST_MODULE "test_parse_range"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/utility.hpp>
#include <cstdint>

BOOST_AUTO_TEST_CASE(test_parse_range_close_close)
{
    {
        const auto actual = mjolnir::parse_range("[1, 10]");
        const auto expect = std::vector<std::size_t>{1,2,3,4,5,6,7,8,9,10};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
    {
        const auto actual = mjolnir::parse_range("[ 1, 10 ]");
        const auto expect = std::vector<std::size_t>{1,2,3,4,5,6,7,8,9,10};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
    {
        const auto actual = mjolnir::parse_range("[1,10]");
        const auto expect = std::vector<std::size_t>{1,2,3,4,5,6,7,8,9,10};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_parse_range_close_open)
{
    {
        const auto actual = mjolnir::parse_range("[1, 10)");
        const auto expect = std::vector<std::size_t>{1,2,3,4,5,6,7,8,9};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
    {
        const auto actual = mjolnir::parse_range("[ 1, 10 )");
        const auto expect = std::vector<std::size_t>{1,2,3,4,5,6,7,8,9};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
    {
        const auto actual = mjolnir::parse_range("[1,10)");
        const auto expect = std::vector<std::size_t>{1,2,3,4,5,6,7,8,9};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_parse_range_open_close)
{
    {
        const auto actual = mjolnir::parse_range("(1, 10]");
        const auto expect = std::vector<std::size_t>{2,3,4,5,6,7,8,9,10};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
    {
        const auto actual = mjolnir::parse_range("( 1, 10 ]");
        const auto expect = std::vector<std::size_t>{2,3,4,5,6,7,8,9,10};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
    {
        const auto actual = mjolnir::parse_range("(1,10]");
        const auto expect = std::vector<std::size_t>{2,3,4,5,6,7,8,9,10};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
}
BOOST_AUTO_TEST_CASE(test_parse_range_open_open)
{
    {
        const auto actual = mjolnir::parse_range("(1, 10)");
        const auto expect = std::vector<std::size_t>{2,3,4,5,6,7,8,9};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
    {
        const auto actual = mjolnir::parse_range("( 1, 10 )");
        const auto expect = std::vector<std::size_t>{2,3,4,5,6,7,8,9};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
    {
        const auto actual = mjolnir::parse_range("(1,10)");
        const auto expect = std::vector<std::size_t>{2,3,4,5,6,7,8,9};
        BOOST_TEST(expect.size() == actual.size());
        if(expect.size() == actual.size())
        {
            BOOST_TEST(std::equal(expect.begin(), expect.end(), actual.begin()));
        }
    }
}
