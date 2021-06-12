#define BOOST_TEST_MODULE "test_unlimited_cell_list"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/stub_potential.hpp>

#include <mjolnir/util/empty.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/range.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Topology.hpp>
#include <numeric>
#include <random>

BOOST_AUTO_TEST_CASE(test_CellList_UnlimitedBoundary)
{
    namespace test = mjolnir::test;
    mjolnir::LoggerManager::set_default_logger("test_cell_list.log");

    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = test::StubPotential<real_type>;

    constexpr std::size_t N = 1000;
    constexpr double      L = 10.0;
    constexpr double cutoff = 2.0;
    constexpr double margin = 0.25; // threshold = 2.0 * (1+0.25) = 2.5
    constexpr double threshold = cutoff * (1.0 + margin);

    const auto distribute_particle = [](std::mt19937& mt, double l) -> coordinate_type
    {
        return coordinate_type(
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt)
        );
    };

    std::vector<std::size_t> participants(N);
    std::iota(participants.begin(), participants.end(), 0u);

    test::StubParameterList<traits_type> params(cutoff, participants, {},
            typename test::StubParameterList<traits_type>::ignore_molecule_type("Nothing"),
            typename test::StubParameterList<traits_type>::ignore_group_type   ({})
            );

    mjolnir::System<traits_type> sys(N, boundary_type{});
    mjolnir::Topology topol(N);

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    topol.construct_molecules();


    mjolnir::SpatialPartition<traits_type, potential_type> vlist(mjolnir::make_unique<
        mjolnir::UnlimitedGridCellList<traits_type, potential_type>>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, params);
    vlist.make(sys, params);
    BOOST_TEST(vlist.valid());

    for(const auto i : params.leading_participants())
    {
        for(std::size_t j=i+1; j<N; ++j)
        {
            const auto partners = vlist.partners(i);
            if(std::find_if(partners.begin(), partners.end(),
                [=](const neighbor_type& elem) -> bool {return elem.index == j;}
                        ) == partners.end())
            {
                // should be enough distant (>= threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(i), sys.position(j)));
                BOOST_TEST(dist >= threshold);
            }
            else
            {
                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(i), sys.position(j)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_CellList_UnlimitedBoundary_partial)
{
    namespace test = mjolnir::test;
    mjolnir::LoggerManager::set_default_logger("test_cell_list.log");

    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = test::StubPotential<real_type>;

    constexpr std::size_t N = 1000;
    constexpr double      L = 10.0;
    constexpr double cutoff = 2.0;
    constexpr double margin = 0.25; // threshold = 2.0 * (1+0.25) = 2.5
    constexpr double threshold = cutoff * (1.0 + margin);

    const auto distribute_particle = [](std::mt19937& mt, double l) -> coordinate_type
    {
        return coordinate_type(
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt)
        );
    };

    std::vector<std::size_t> participants(500);
    // [200 ... 699]
    std::iota(participants.begin(), participants.end(), 200u);

    test::StubParameterList<traits_type> params(cutoff, participants, {},
            typename test::StubParameterList<traits_type>::ignore_molecule_type("Nothing"),
            typename test::StubParameterList<traits_type>::ignore_group_type   ({})
            );

    mjolnir::System<traits_type> sys(N, boundary_type{});
    mjolnir::Topology topol(N);

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    topol.construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(mjolnir::make_unique<
        mjolnir::UnlimitedGridCellList<traits_type, potential_type>>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, params);
    vlist.make(sys, params);
    BOOST_TEST(vlist.valid());

    for(const auto i : params.leading_participants())
    {
        const auto partners = vlist.partners(i);

        // if particle i is not related to the potential, it should not have
        // any interacting partners.
        if(params.participants().end() == std::find(
            params.participants().begin(), params.participants().end(), i))
        {
            BOOST_TEST(partners.size() == 0u);
            continue;
        }

        for(std::size_t j=i+1; j<N; ++j)
        {
            if(std::find_if(partners.begin(), partners.end(),
                [=](const neighbor_type& elem) -> bool {return elem.index == j;}
                        ) == partners.end())
            {
                // should be enough distant (>= threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(i), sys.position(j)));
                const bool enough_distant = dist >= threshold;

                // or not a participant
                const bool is_participant = params.participants().end() != std::find(
                    params.participants().begin(), params.participants().end(), j);

                const bool is_ok = enough_distant || (!is_participant);
                BOOST_TEST(is_ok);
            }
            else
            {
                // should be a participant
                const bool found = params.participants().end() != std::find(
                    params.participants().begin(), params.participants().end(), j);
                BOOST_TEST(found);

                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(i), sys.position(j)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_CellList_UnlimitedBoundary_partial_2)
{
    namespace test = mjolnir::test;
    mjolnir::LoggerManager::set_default_logger("test_cell_list.log");

    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = test::StubPotential<real_type>;

    constexpr std::size_t N = 1000;
    constexpr double      L = 10.0;
    constexpr double cutoff = 2.0;
    constexpr double margin = 0.25; // threshold = 2.0 * (1+0.25) = 2.5
    constexpr double threshold = cutoff * (1.0 + margin);

    const auto distribute_particle = [](std::mt19937& mt, double l) -> coordinate_type
    {
        return coordinate_type(
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt),
            l * std::generate_canonical<real_type, std::numeric_limits<real_type>::digits>(mt)
        );
    };

    std::vector<std::size_t> participants; participants.reserve(500);
    for(std::size_t i=0; i<500; ++i)
    {
        participants.push_back(i*2);
    }

    test::StubParameterList<traits_type> params(cutoff, participants, {},
            typename test::StubParameterList<traits_type>::ignore_molecule_type("Nothing"),
            typename test::StubParameterList<traits_type>::ignore_group_type   ({})
            );

    mjolnir::System<traits_type> sys(N, boundary_type{});
    mjolnir::Topology topol(N);

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    topol.construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(mjolnir::make_unique<
        mjolnir::UnlimitedGridCellList<traits_type, potential_type>>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, params);
    vlist.make(sys, params);
    BOOST_TEST(vlist.valid());

    for(const auto i : params.leading_participants())
    {
        const auto partners = vlist.partners(i);

        // if particle i is not related to the potential, it should not have
        // any interacting partners.
        if(params.participants().end() == std::find(
            params.participants().begin(), params.participants().end(), i))
        {
            BOOST_TEST(partners.size() == 0u);
            continue;
        }

        for(std::size_t j=i+1; j<N; ++j)
        {
            if(std::find_if(partners.begin(), partners.end(),
                [=](const neighbor_type& elem) -> bool {return elem.index == j;}
                        ) == partners.end())
            {
                // should be enough distant (>= threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(i), sys.position(j)));
                const bool enough_distant = dist >= threshold;

                // or not a participant
                const bool is_participant = params.participants().end() != std::find(
                    params.participants().begin(), params.participants().end(), j);

                const bool is_ok = enough_distant || (!is_participant);
                BOOST_TEST(is_ok);
            }
            else
            {
                // should be a participant
                const bool found = params.participants().end() != std::find(
                    params.participants().begin(), params.participants().end(), j);
                BOOST_TEST(found);

                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(i), sys.position(j)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_CellList_UnlimitedBoundary_clone)
{
    mjolnir::LoggerManager::set_default_logger("test_cell_list.log");

    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = typename traits_type::real_type;
    using potential_type  = test::StubPotential<real_type>;

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(mjolnir::make_unique<
        mjolnir::UnlimitedGridCellList<traits_type, potential_type>>(10.0));

    mjolnir::SpatialPartition<traits_type, potential_type> vlist2(vlist);

    BOOST_TEST(vlist.margin() == vlist2.margin());
}
