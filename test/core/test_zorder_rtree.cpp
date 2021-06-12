#define BOOST_TEST_MODULE "test_zorder_rtree"

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
#include <mjolnir/core/ZorderRTree.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Topology.hpp>
#include <random>

using traits_to_be_tested = std::tuple<
    mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>,
    mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>,
    mjolnir::SimulatorTraits<float, mjolnir::UnlimitedBoundary>,
    mjolnir::SimulatorTraits<float, mjolnir::CuboidalPeriodicBoundary>
>;

template<typename Boundary, typename Coordinate>
typename std::enable_if<mjolnir::is_unlimited_boundary<Boundary>::value, Boundary>::type
get_boundary(const Coordinate&, const Coordinate&) noexcept
{
    return Boundary{};
}
template<typename Boundary, typename Coordinate>
typename std::enable_if<mjolnir::is_cuboidal_periodic_boundary<Boundary>::value, Boundary>::type
get_boundary(const Coordinate& lower, const Coordinate& higher) noexcept
{
    return Boundary{lower, higher};
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_ZorderRTree_full_interaction, traits_type, traits_to_be_tested)
{
    mjolnir::LoggerManager::set_default_logger("test_zorder_rtree.log");
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = test::StubPotential<real_type>;

    constexpr std::size_t N = 1000;
    constexpr double      L = 10.0;
    constexpr double cutoff = 1.0;
    constexpr double margin = 0.5;
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

    mjolnir::System<traits_type> sys(N, get_boundary<boundary_type>(
                coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));
    mjolnir::Topology topol(N);

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    topol.construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(
        mjolnir::make_unique<mjolnir::ZorderRTree<traits_type, potential_type>>(margin));
    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : pot.leading_participants())
    {
        for(std::size_t j=i+1; j<N; ++j)
        {
            const auto partners = vlist.partners(i);
            if(std::find_if(partners.begin(), partners.end(),
                    [=](const neighbor_type& elem){return elem.index == j;}
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

BOOST_AUTO_TEST_CASE_TEMPLATE(test_ZorderRTree_partial_interaction, traits_type, traits_to_be_tested)
{
    mjolnir::LoggerManager::set_default_logger("test_zorder_rtree.log");
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

    mjolnir::System<traits_type> sys(N, get_boundary<boundary_type>(
                coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));
    mjolnir::Topology topol(N);

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    topol.construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(
        mjolnir::make_unique<mjolnir::ZorderRTree<traits_type, potential_type>>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : pot.leading_participants())
    {
        const auto partners = vlist.partners(i);

        // if particle i is not related to the potential, it should not have
        // any interacting partners.
        if(pot.participants().end() == std::find(
            pot.participants().begin(), pot.participants().end(), i))
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
                const bool is_participant = pot.participants().end() != std::find(
                    pot.participants().begin(), pot.participants().end(), j);

                const bool is_ok = enough_distant || (!is_participant);
                BOOST_TEST(is_ok);
            }
            else
            {
                // should be a participant
                const bool found = pot.participants().end() != std::find(
                    pot.participants().begin(), pot.participants().end(), j);
                BOOST_TEST(found);

                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(i), sys.position(j)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_ZorderRTree_partial_interaction_incontiguous, traits_type, traits_to_be_tested)
{
    mjolnir::LoggerManager::set_default_logger("test_zorder_rtree.log");
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
        participants.push_back(i * 2);
    }

    test::StubParameterList<traits_type> params(cutoff, participants, {},
            typename test::StubParameterList<traits_type>::ignore_molecule_type("Nothing"),
            typename test::StubParameterList<traits_type>::ignore_group_type   ({})
            );


    mjolnir::System<traits_type> sys(N, get_boundary<boundary_type>(
                coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));
    mjolnir::Topology topol(N);

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    topol.construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(
            mjolnir::make_unique<mjolnir::ZorderRTree<traits_type, potential_type>>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : pot.leading_participants())
    {
        const auto partners = vlist.partners(i);

        // if particle i is not related to the potential, it should not have
        // any interacting partners.
        if(pot.participants().end() == std::find(
            pot.participants().begin(), pot.participants().end(), i))
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
                const bool is_participant = pot.participants().end() != std::find(
                    pot.participants().begin(), pot.participants().end(), j);

                const bool is_ok = enough_distant || (!is_participant);
                BOOST_TEST(is_ok);
            }
            else
            {
                // should be a participant
                const bool found = pot.participants().end() != std::find(
                    pot.participants().begin(), pot.participants().end(), j);
                BOOST_TEST(found);

                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(i), sys.position(j)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_ZorderRTree_clone, traits_type, traits_to_be_tested)
{
    mjolnir::LoggerManager::set_default_logger("test_zorder_rtree.log");
    using real_type       = typename traits_type::real_type;
    using potential_type  = test::StubPotential<real_type>;

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(
        mjolnir::make_unique<mjolnir::ZorderRTree<traits_type, potential_type>>(1.0));

    mjolnir::SpatialPartition<traits_type, potential_type> vlist2(vlist);

    BOOST_TEST(vlist.margin() == vlist2.margin());
}
