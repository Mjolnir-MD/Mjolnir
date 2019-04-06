#define BOOST_TEST_MODULE "test_range"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/range.hpp>

BOOST_AUTO_TEST_CASE(test_range)
{
    {
        const std::vector<int> is(10, 42);
        const auto rg = mjolnir::make_range(is.begin(), is.end());

        BOOST_TEST(rg.size() == 10u);

        for(const int i : rg)
        {
            BOOST_TEST(i == 42);
        }
    }

    {
        const std::vector<int> is = {0,1,2,3,4,5,6,7,8,9};
        const auto rg = mjolnir::make_range(is.begin(), is.end());

        int a = 0;
        for(const int i : rg)
        {
            BOOST_TEST(i == a);
            ++a;
        }
    }
}
