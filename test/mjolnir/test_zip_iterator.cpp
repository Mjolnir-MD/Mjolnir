#define BOOST_TEST_MODULE "test_zip_iterator"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <mjolnir/util/range.hpp>

BOOST_AUTO_TEST_CASE(test_zipped_traverse)
{
    {
        const std::vector<int>    is(10, 42);
        const std::vector<double> ds(10, 3.14);

        for(const auto& i : mjolnir::make_range(
                mjolnir::make_zip(is.begin(), ds.begin()),
                mjolnir::make_zip(is.end(),   ds.end())))
        {
            BOOST_CHECK_EQUAL(std::get<0>(i), 42);
            BOOST_CHECK_EQUAL(std::get<1>(i), 3.14);
        }
    }

    {
        const std::vector<int> is = {0,1,2,3,4,5,6,7,8,9};
        const std::vector<int> rv = {9,8,7,6,5,4,3,2,1,0};

        std::size_t a = 0;
        for(const auto& i : mjolnir::make_range(
                mjolnir::make_zip(is.begin(), rv.begin()),
                mjolnir::make_zip(is.end(),   rv.end())))
        {
            BOOST_CHECK_EQUAL(std::get<0>(i), a);
            BOOST_CHECK_EQUAL(std::get<1>(i), 9 - a);
            ++a;
        }
    }
}

BOOST_AUTO_TEST_CASE(test_zipped_substitution)
{
    {
        const std::vector<int> ref1 = {0,1,2,3,4,5,6,7,8,9};
        const std::vector<int> ref2 = {9,8,7,6,5,4,3,2,1,0};

        std::vector<int> ans1(ref1.size());
        std::vector<int> ans2(ref2.size());

        std::copy(mjolnir::make_zip(ref1.begin(), ref2.begin()),
                  mjolnir::make_zip(ref1.end(),   ref2.end()),
                  mjolnir::make_zip(ans1.begin(), ans2.begin()));

        for(std::size_t i=0; i<ref1.size(); ++i)
        {
            BOOST_CHECK_EQUAL(ref1.at(i), ans1.at(i));
            BOOST_CHECK_EQUAL(ref2.at(i), ans2.at(i));
        }
    }
}
