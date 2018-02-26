#define BOOST_TEST_MODULE "test_range"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/range.hpp>

BOOST_AUTO_TEST_CASE(test_range)
{
    {
        const std::vector<int> is(10, 42);
        const auto rg = mjolnir::make_range(is.begin(), is.end());

        BOOST_CHECK_EQUAL(rg.size(), 10);

        for(const auto& i : rg)
        {
            BOOST_CHECK_EQUAL(i, 42);
        }
    }

    {
        const std::vector<int> is = {0,1,2,3,4,5,6,7,8,9};
        const auto rg = mjolnir::make_range(is.begin(), is.end());

        std::size_t a = 0;
        for(const auto& i : rg)
        {
            BOOST_CHECK_EQUAL(i, a);
            ++a;
        }
    }
}
