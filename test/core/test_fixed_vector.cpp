#define BOOST_TEST_MODULE "test_fixed_vector"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/fixed_vector.hpp>

BOOST_AUTO_TEST_CASE(test_initialization)
{
    {
        mjolnir::fixed_vector<int, 8> fv;
        BOOST_TEST(fv.empty());
        BOOST_TEST(fv.size() == 0u);
        BOOST_TEST(static_cast<bool>(fv.begin()   == fv.end()  ));
        BOOST_TEST(static_cast<bool>(fv.cbegin()  == fv.cend() ));
        BOOST_TEST(static_cast<bool>(fv.rbegin()  == fv.rend() ));
        BOOST_TEST(static_cast<bool>(fv.crbegin() == fv.crend()));
    }

    {
        mjolnir::fixed_vector<int, 8> fv{1, 2, 3};
        BOOST_TEST( ! fv.empty());
        BOOST_TEST(fv.size()  == 3u);
        BOOST_TEST(fv.front() == 1);
        BOOST_TEST(fv.back()  == 3);

        BOOST_TEST(fv.at(0)  == 1);
        BOOST_TEST(fv.at(1)  == 2);
        BOOST_TEST(fv.at(2)  == 3);
        BOOST_TEST(fv[0]     == 1);
        BOOST_TEST(fv[1]     == 2);
        BOOST_TEST(fv[2]     == 3);

        int i=1;
        for(const auto elem : fv)
        {
            BOOST_TEST(i == elem);
            i++;
        }
    }

    {
        mjolnir::fixed_vector<int, 8> fv;
        fv.push_back(1);
        fv.push_back(2);
        fv.push_back(3);

        BOOST_TEST( ! fv.empty());
        BOOST_TEST(fv.size()  == 3u);
        BOOST_TEST(fv.front() == 1);
        BOOST_TEST(fv.back()  == 3);

        BOOST_TEST(fv.at(0)  == 1);
        BOOST_TEST(fv.at(1)  == 2);
        BOOST_TEST(fv.at(2)  == 3);
        BOOST_TEST(fv[0]     == 1);
        BOOST_TEST(fv[1]     == 2);
        BOOST_TEST(fv[2]     == 3);

        int i=1;
        for(const auto elem : fv)
        {
            BOOST_TEST(i == elem);
            i++;
        }
    }
}
