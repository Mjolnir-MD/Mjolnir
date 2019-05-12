#define BOOST_TEST_MODULE "test_omp_sort"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/omp/sort.hpp>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>


BOOST_AUTO_TEST_CASE(test_omp_sort_large)
{
    BOOST_TEST_WARN(omp_get_max_threads() > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    const std::size_t N = 10000;
    for(std::size_t i=0, e=omp_get_max_threads(); i<e; ++i)
    {
        std::vector<int> vec(N+i, 0);
        std::vector<int> buf(N+i, 0);

        std::iota(vec.begin(), vec.end(), 1);
        std::mt19937 mt(123456789);
        std::shuffle(vec.begin(), vec.end(), mt);

        mjolnir::omp::sort(vec, buf);
        BOOST_TEST(std::is_sorted(vec.begin(), vec.end()));
    }
}

BOOST_AUTO_TEST_CASE(test_omp_sort_cmp_large)
{
    BOOST_TEST_WARN(omp_get_max_threads() > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    const std::size_t N = 10000;
    const auto comp = [](const std::pair<int, int>& lhs,
                         const std::pair<int, int>& rhs) -> bool {
        return lhs.first < rhs.first;
    };

    for(std::size_t i=0, e=omp_get_max_threads(); i<e; ++i)
    {
        std::vector<std::pair<int, int>> vec(N+i, std::make_pair(0, 0));
        std::vector<std::pair<int, int>> buf(N+i, std::make_pair(0, 0));
        for(std::size_t j=0; j<vec.size(); ++j)
        {
            vec.at(j).first = j+1;
        }

        std::mt19937 mt(123456789);
        std::shuffle(vec.begin(), vec.end(), mt);

        mjolnir::omp::sort(vec, buf, comp);
        BOOST_TEST(std::is_sorted(vec.begin(), vec.end(), comp));
    }
}

// test if a number of elements in a vector is less than the number of threads
BOOST_AUTO_TEST_CASE(test_omp_sort_small)
{
    BOOST_TEST_WARN(omp_get_max_threads() > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    const std::size_t N = 2;
    for(std::size_t i=0, e=omp_get_max_threads(); i<e; ++i)
    {
        std::vector<int> vec(N+i, 0);
        std::vector<int> buf(N+i, 0);

        std::iota(vec.begin(), vec.end(), 1);
        std::mt19937 mt(123456789);
        std::shuffle(vec.begin(), vec.end(), mt);

        mjolnir::omp::sort(vec, buf);
        BOOST_TEST(std::is_sorted(vec.begin(), vec.end()));
    }
}

BOOST_AUTO_TEST_CASE(test_omp_sort_cmp_small)
{
    BOOST_TEST_WARN(omp_get_max_threads() > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    const std::size_t N = 2;
    const auto comp = [](const std::pair<int, int>& lhs,
                         const std::pair<int, int>& rhs) -> bool {
        return lhs.first < rhs.first;
    };

    for(std::size_t i=0, e=omp_get_max_threads(); i<e; ++i)
    {
        std::vector<std::pair<int, int>> vec(N+i, std::make_pair(0, 0));
        std::vector<std::pair<int, int>> buf(N+i, std::make_pair(0, 0));
        for(std::size_t j=0; j<vec.size(); ++j)
        {
            vec.at(j).first = j+1;
        }

        std::mt19937 mt(123456789);
        std::shuffle(vec.begin(), vec.end(), mt);

        mjolnir::omp::sort(vec, buf, comp);
        BOOST_TEST(std::is_sorted(vec.begin(), vec.end(), comp));
    }
}
