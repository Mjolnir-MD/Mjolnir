#define BOOST_TEST_MODULE "test_neighbor_list"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/empty.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/NeighborList.hpp>

BOOST_AUTO_TEST_CASE(test_NeighborList_no_parameter)
{
    mjolnir::LoggerManager::set_default_logger("test_neighbor_list.log");
    using parameter_type  = mjolnir::empty_t; // XXX

    mjolnir::NeighborList<parameter_type> nlist;
    using neighbor_type =
        typename mjolnir::NeighborList<parameter_type>::neighbor_type;

    static_assert(sizeof(neighbor_type) == sizeof(std::size_t), "");

    // construct dummy neighbor list
    const std::size_t N = 100;
    const std::size_t M =  10;
    for(std::size_t i=0; i<N; ++i)
    {
        std::vector<neighbor_type> partner;
        for(std::size_t j = i + 1; j < i + 1 + M; ++j)
        {
            partner.emplace_back(j, parameter_type{});
        }
        nlist.add_list_for(i, partner.begin(), partner.end());
    }

    // check the list is correctly built
    for(std::size_t i=0; i<N; ++i)
    {
        const auto partners = nlist.at(i);
        BOOST_TEST(partners.size() == M);
        for(std::size_t k=0; k<M; ++k)
        {
            BOOST_TEST(partners.at(k).index == i + 1 + k);
            // no partner here.
        }
    }
}

BOOST_AUTO_TEST_CASE(test_NeighborList_with_parameter)
{
    mjolnir::LoggerManager::set_default_logger("test_neighbor_list.log");
    // contain this parameter to neighbor list
    using parameter_type  = std::pair<std::size_t, std::size_t>;

    mjolnir::NeighborList<parameter_type> nlist;
    using neighbor_type =
        typename mjolnir::NeighborList<parameter_type>::neighbor_type;

    static_assert(
        sizeof(neighbor_type) == sizeof(std::size_t) + sizeof(parameter_type),
        "");

    // construct dummy neighbor list
    const std::size_t N = 100;
    const std::size_t M =  10;
    for(std::size_t i=0; i<N; ++i)
    {
        std::vector<neighbor_type> partner;
        for(std::size_t j = i + 1; j < i + 1 + M; ++j)
        {
            partner.emplace_back(j, parameter_type{i, j});
        }
        nlist.add_list_for(i, partner.begin(), partner.end());
    }

    // check the list is correctly built
    for(std::size_t i=0; i<N; ++i)
    {
        const auto partners = nlist.at(i);
        BOOST_TEST(partners.size() == M);
        for(std::size_t k=0; k<M; ++k)
        {
            BOOST_TEST(partners.at(k).index              == i + 1 + k);
            BOOST_TEST(partners.at(k).potential().first  == i        );
            BOOST_TEST(partners.at(k).potential().second == i + 1 + k);
        }
    }
}
