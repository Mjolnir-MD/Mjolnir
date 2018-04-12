#define BOOST_TEST_MODULE "test_split_string"
#include <boost/test/included/unit_test.hpp>
#include <mjolnir/util/string.hpp>
#include <jarngreipr/util/split_string.hpp>

using mjolnir::operator"" _str;

const auto source = "ABC-DE-FGHI-J-K-L"_str;
const std::array<std::string, 6> answer = {{
    "ABC"_str, "DE"_str, "FGHI"_str, "J"_str, "K"_str, "L"_str
}};

BOOST_AUTO_TEST_CASE(test_split_string_range_based_for)
{
    std::size_t idx = 0;
    for(auto&& elem : mjolnir::split_string(source, '-'))
    {
        BOOST_CHECK(idx < 6);
        BOOST_CHECK_EQUAL(elem, answer.at(idx));
        ++idx;
    }
}

BOOST_AUTO_TEST_CASE(test_split_string_increment_decrement)
{
    const auto range = mjolnir::split_string(source, '-');
    auto iter = range.begin();
    BOOST_CHECK_EQUAL(*iter, answer.at(0)); ++iter;
    BOOST_CHECK_EQUAL(*iter, answer.at(1)); ++iter;
    BOOST_CHECK_EQUAL(*iter, answer.at(2)); ++iter;
    BOOST_CHECK_EQUAL(*iter, answer.at(3)); ++iter;
    BOOST_CHECK_EQUAL(*iter, answer.at(4)); ++iter;
    BOOST_CHECK_EQUAL(*iter, answer.at(5));
    --iter; BOOST_CHECK_EQUAL(*iter, answer.at(4));
    --iter; BOOST_CHECK_EQUAL(*iter, answer.at(3));
    --iter; BOOST_CHECK_EQUAL(*iter, answer.at(2));
    --iter; BOOST_CHECK_EQUAL(*iter, answer.at(1));
    --iter; BOOST_CHECK_EQUAL(*iter, answer.at(0));
    BOOST_CHECK(iter == range.begin());
}

