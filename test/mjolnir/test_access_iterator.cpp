#define BOOST_TEST_MODULE "test_access_iterator"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/util/access_iterator.hpp>
#include <mjolnir/util/range.hpp>

BOOST_AUTO_TEST_CASE(test_immutable_access_iterator)
{
    const auto a = [](const std::pair<int, double>& p) noexcept -> const int& {
        return p.first;
    };

    const std::vector<std::pair<int, double>> vec ={
        {1,  3.14},
        {2,  6.28},
        {3,  9.42},
        {4, 12.56},
    };
    using iter = mjolnir::access_iterator<decltype(vec.cbegin()), decltype(a)>;

    std::size_t idx=0;
    for(iter i(vec.cbegin(), a), e(vec.cend(), a); i!=e; ++i)
    {
        BOOST_CHECK_EQUAL(*i, vec.at(idx).first);
        BOOST_CHECK_EQUAL(*i, a(vec.at(idx)));
        ++idx;
    }
}

BOOST_AUTO_TEST_CASE(test_mutable_access_iterator)
{
    const auto a = [](std::pair<int, double>& p) noexcept -> int& {
        return p.first;
    };

    std::vector<std::pair<int, double>> vec ={
        {1,  3.14},
        {2,  6.28},
        {3,  9.42},
        {4, 12.56},
    };
    using iter = mjolnir::access_iterator<decltype(vec.begin()), decltype(a)>;

    std::size_t idx=0;
    for(iter i(vec.begin(), a), e(vec.end(), a); i!=e; ++i)
    {
        BOOST_CHECK_EQUAL(*i, vec.at(idx).first);
        BOOST_CHECK_EQUAL(*i, a(vec.at(idx)));
        ++idx;
    }

    for(iter i(vec.begin(), a), e(vec.end(), a); i!=e; ++i)
    {
        *i = 5;
    }

    idx = 0;
    for(iter i(vec.begin(), a), e(vec.end(), a); i!=e; ++i)
    {
        BOOST_CHECK_EQUAL(*i, 5);
        BOOST_CHECK_EQUAL(*i, a(vec.at(idx)));
        ++idx;
    }
}

BOOST_AUTO_TEST_CASE(test_make_access_iterator)
{
    const auto a = [](const std::pair<int, double>& p) noexcept -> const int& {
        return p.first;
    };

    const std::vector<std::pair<int, double>> vec ={
        {1,  3.14},
        {2,  6.28},
        {3,  9.42},
        {4, 12.56},
    };

    std::size_t idx=0;
    for(auto i(mjolnir::make_access_iterator(vec.cbegin(), a)),
             e(mjolnir::make_access_iterator(vec.cend(), a)); i!=e; ++i)
    {
        BOOST_CHECK_EQUAL(*i, vec.at(idx).first);
        BOOST_CHECK_EQUAL(*i, a(vec.at(idx)));
        ++idx;
    }
}

BOOST_AUTO_TEST_CASE(test_make_access_iterator_and_range)
{
    const auto a = [](const std::pair<int, double>& p) noexcept -> const int& {
        return p.first;
    };

    const std::vector<std::pair<int, double>> vec ={
        {1,  3.14},
        {2,  6.28},
        {3,  9.42},
        {4, 12.56},
    };

    const auto range = mjolnir::make_range(
            mjolnir::make_access_iterator(vec.cbegin(), a),
            mjolnir::make_access_iterator(vec.cend(), a));
    BOOST_CHECK_EQUAL(range.size(), 4u);

    std::size_t idx=0;
    for(auto i : range)
    {
        BOOST_CHECK_EQUAL(i, vec.at(idx).first);
        BOOST_CHECK_EQUAL(i, a(vec.at(idx)));
        ++idx;
    }
}
