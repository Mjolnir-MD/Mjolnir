#define BOOST_TEST_MODULE "test_split_string"
#include <boost/test/included/unit_test.hpp>
#include <mjolnir/util/string.hpp>
#include <jarngreipr/io/read_chain_ids.hpp>

BOOST_AUTO_TEST_CASE(test_split_string)
{
    using mjolnir::operator"" _str;
    {
        const auto splitted = jarngreipr::split_string("ABC-DEF-GHIJ"_str, '-');
        BOOST_CHECK_EQUAL(splitted.size(), 3);
        BOOST_CHECK_EQUAL(splitted.at(0), "ABC"_str);
        BOOST_CHECK_EQUAL(splitted.at(1), "DEF"_str);
        BOOST_CHECK_EQUAL(splitted.at(2), "GHIJ"_str);
    }

    {
        const auto splitted = jarngreipr::split_string("A-B-C"_str, '-');
        BOOST_CHECK_EQUAL(splitted.size(), 3);
        BOOST_CHECK_EQUAL(splitted.at(0), "A"_str);
        BOOST_CHECK_EQUAL(splitted.at(1), "B"_str);
        BOOST_CHECK_EQUAL(splitted.at(2), "C"_str);
    }
}

BOOST_AUTO_TEST_CASE(test_read_chain_ids)
{
    using mjolnir::operator"" _str;

    {
        const auto splitted = jarngreipr::read_chain_ids("D"_str);
        BOOST_CHECK_EQUAL(splitted.size(), 1);
        BOOST_CHECK_EQUAL(splitted.at(0), 'D');
    }

    {
        const auto splitted = jarngreipr::read_chain_ids("B-C"_str);
        BOOST_CHECK_EQUAL(splitted.size(), 2);
        BOOST_CHECK_EQUAL(splitted.at(0), 'B');
        BOOST_CHECK_EQUAL(splitted.at(1), 'C');
    }

    {
        const auto splitted = jarngreipr::read_chain_ids("B-D"_str);
        BOOST_CHECK_EQUAL(splitted.size(), 3);
        BOOST_CHECK_EQUAL(splitted.at(0), 'B');
        BOOST_CHECK_EQUAL(splitted.at(1), 'C');
        BOOST_CHECK_EQUAL(splitted.at(2), 'D');
    }

    {
        const auto splitted = jarngreipr::read_chain_ids("B&D"_str);
        BOOST_CHECK_EQUAL(splitted.size(), 2);
        BOOST_CHECK_EQUAL(splitted.at(0), 'B');
        BOOST_CHECK_EQUAL(splitted.at(1), 'D');
    }

    {
        const auto splitted = jarngreipr::read_chain_ids("A-C&F-G"_str);
        BOOST_CHECK_EQUAL(splitted.size(), 5);
        BOOST_CHECK_EQUAL(splitted.at(0), 'A');
        BOOST_CHECK_EQUAL(splitted.at(1), 'B');
        BOOST_CHECK_EQUAL(splitted.at(2), 'C');
        BOOST_CHECK_EQUAL(splitted.at(3), 'F');
        BOOST_CHECK_EQUAL(splitted.at(4), 'G');
    }
}
