#define BOOST_TEST_MODULE "test_topology"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/Topology.hpp>
#include <cstdint>

BOOST_AUTO_TEST_CASE(topology_list_adjacent_is_unique)
{
    std::size_t N = 100;
    mjolnir::Topology top(N);
    for(std::size_t i=0; i<N-1; ++i)
    {
        top.add_connection(i, i+1, "bond");
    }

    for(std::size_t i=0; i<N; ++i)
    {
        auto tmp = top.list_adjacent_within(i, 3, "bond");
        std::sort(tmp.begin(), tmp.end());
        const auto result = std::unique(tmp.begin(), tmp.end());
        BOOST_CHECK(result == tmp.end());
    }
}

BOOST_AUTO_TEST_CASE(topology_list_within)
{
    std::size_t N = 100;
    mjolnir::Topology top(100);
    for(std::size_t i=0; i<N-1; ++i)
    {
        top.add_connection(i, i+1, "bond");
    }

    for(std::size_t i=0; i<N; ++i)
    {
        const auto adj3 = top.list_adjacent_within(i, 3, "bond");
        for(const auto t : adj3)
        {
            BOOST_CHECK(0 <= t && t < N);
            const auto i_ = static_cast<std::int32_t>(i);
            const auto t_ = static_cast<std::int32_t>(t);
            const auto diff = std::abs(i_ - t_);
            BOOST_CHECK(diff <= 3);
        }

        if(4 <= i && i <= N-4)
        {
            BOOST_CHECK_EQUAL(adj3.size(), 7);
        }

        const auto adj_b3 = top.list_adjacent_within(i, 3, "bond");
        BOOST_CHECK_EQUAL(adj3.size(), adj_b3.size());

        if(adj3.size() == adj_b3.size())
        {
            const auto adj3_is_equal_to_adj_b3 = std::is_permutation(
                    adj3.begin(), adj3.end(), adj_b3.begin());
            BOOST_CHECK(adj3_is_equal_to_adj_b3);
        }

        const auto adj_c3 = top.list_adjacent_within(i, 3, "contact");
        BOOST_CHECK_EQUAL(adj_c3.size(), 1);
        if(!adj_c3.empty())
        {
            BOOST_CHECK_EQUAL(adj_c3.front(), i);
        }
    }
}

BOOST_AUTO_TEST_CASE(topology_contacts_list_within)
{
    std::size_t N = 100;
    mjolnir::Topology top(100);
    for(std::size_t i=0; i<N-1; ++i)
    {
        top.add_connection(i, i+1, "bond");
    }
    for(std::size_t i=0; i<N-4; ++i)
    {
        top.add_connection(i, i+4, "contact");
    }

    for(std::size_t i=0; i<N; ++i)
    {
        const auto adj3 = top.list_adjacent_within(i, 3, "bond");
        for(const auto t : adj3)
        {
            BOOST_CHECK(0 <= t && t < N);
            const auto i_ = static_cast<std::int32_t>(i);
            const auto t_ = static_cast<std::int32_t>(t);
            const auto diff = std::abs(i_ - t_);
            BOOST_CHECK(diff <= 3);
        }

        if(4 <= i && i <= N-4)
        {
            BOOST_CHECK_EQUAL(adj3.size(), 7);
        }
    }
    for(std::size_t i=0; i<N; ++i)
    {
        const auto bd = top.list_adjacent_within(i, 1, "bond");

        for(const auto a : bd)
        {
            BOOST_CHECK(a == i ||  (top.has_connection(i, a, "bond")));
            BOOST_CHECK(a == i || !(top.has_connection(i, a, "contact")));
        }

        const auto ct = top.list_adjacent_within(i, 1, "contact");
        for(const auto a : ct)
        {
            BOOST_CHECK(a == i || !(top.has_connection(i, a, "bond")));
            BOOST_CHECK(a == i ||  (top.has_connection(i, a, "contact")));
        }
    }
}


BOOST_AUTO_TEST_CASE(topology_consistency_has_connection_and_list_within)
{
    std::size_t N = 100;
    mjolnir::Topology top(N);
    for(std::size_t i=0; i<N-1; ++i)
    {
        top.add_connection(i, i+1, "bond");
    }

    for(std::size_t i=0; i<N; ++i)
    {
        BOOST_CHECK((top.has_connection(i, i, "bond")));
        for(const auto adj : top.list_adjacent_within(i, 1, "bond"))
        {
            BOOST_CHECK((top.has_connection(i, adj, "bond")));
            BOOST_CHECK((top.has_connection(adj, i, "bond")));
        }
    }
}

BOOST_AUTO_TEST_CASE(topology_erase_connection)
{
    std::size_t N = 100;
    mjolnir::Topology top(N);
    for(std::size_t i=0; i<N-1; ++i)
    {
        top.add_connection(i, i+1, "bond");
    }

    for(std::size_t i=0; i<N-1; ++i)
    {
        BOOST_CHECK((top.has_connection(i, i+1, "bond")));

        auto tmp = top;
        BOOST_CHECK((tmp.has_connection(i, i+1, "bond")));

        tmp.erase_connection(i, i+1, "bond");
        BOOST_CHECK(!(tmp.has_connection(i, i+1, "bond")));
    }
}
