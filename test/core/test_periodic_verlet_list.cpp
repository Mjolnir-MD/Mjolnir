#define BOOST_TEST_MODULE "test_periodic_verlet_list"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/empty.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Topology.hpp>
#include <random>

// has non-empty parameter type.
template<typename T>
struct dummy_potential
{
    using real_type      = T;
    using parameter_type = std::pair<std::size_t, std::size_t>;
    using pair_parameter_type = parameter_type;

    using topology_type        = mjolnir::Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;

    explicit dummy_potential(const real_type cutoff,
                             const std::vector<std::size_t>& participants)
        : cutoff_(cutoff), participants_(participants)
    {}

    real_type max_cutoff_length() const noexcept {return this->cutoff_;}

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        return parameter_type{i, j};
    }

    bool is_ignored_molecule(std::size_t, std::size_t) const {return false;}
    bool is_ignored_group   (std::string, std::string) const {return false;}
    bool has_interaction    (std::size_t i, std::size_t j) const
    {
        return i != j;
    }

    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within() const
    {
        return std::vector<std::pair<connection_kind_type, std::size_t>>{};
    }

    std::vector<std::size_t> const& participants() const noexcept
    {
        return this->participants_;
    }

    std::string name() const {return "dummy potential";}

    real_type cutoff_;
    std::vector<std::size_t> participants_;
};

BOOST_AUTO_TEST_CASE(test_VerletList_PeriodicBoundary)
{
    mjolnir::LoggerManager::set_default_logger("test_verlet_list.log");
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = dummy_potential<real_type>;

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
    std::iota(participants.begin(), participants.end(), 0);

    dummy_potential<real_type> pot(cutoff, participants);

    mjolnir::System<traits_type> sys(N, boundary_type(coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    sys.topology().construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(mjolnir::make_unique<
        mjolnir::VerletList<traits_type, potential_type>>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : participants)
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
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist >= threshold);
            }
            else
            {
                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist < threshold);
            }
        }
    }

    // check parameter_type.
    for(const auto i : participants)
    {
        for(const auto& p_j : vlist.partners(i))
        {
            const std::size_t j = p_j.index;
            BOOST_TEST(p_j.parameter().first  == pot.prepare_params(i, j).first);
            BOOST_TEST(p_j.parameter().second == pot.prepare_params(i, j).second);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_VerletList_PeriodicBoundary_partial)
{
    mjolnir::LoggerManager::set_default_logger("test_verlet_list.log");
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = dummy_potential<real_type>;

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

    dummy_potential<real_type> pot(cutoff, participants);

    mjolnir::System<traits_type> sys(N, boundary_type(coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    sys.topology().construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(mjolnir::make_unique<
        mjolnir::VerletList<traits_type, potential_type>>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : participants)
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
                            sys.position(j) - sys.position(i)));
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
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_VerletList_PeriodicBoundary_partial_2)
{
    mjolnir::LoggerManager::set_default_logger("test_verlet_list.log");
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = dummy_potential<real_type>;

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

    dummy_potential<real_type> pot(cutoff, participants);

    mjolnir::System<traits_type> sys(N, boundary_type(coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    sys.topology().construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(mjolnir::make_unique<
        mjolnir::VerletList<traits_type, potential_type>>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : participants)
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
                            sys.position(j) - sys.position(i)));
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
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// allow both (i,j) and (j,i)

BOOST_AUTO_TEST_CASE(test_VerletList_PeriodicBoundary_interchange)
{
    mjolnir::LoggerManager::set_default_logger("test_verlet_list.log");
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = dummy_potential<real_type>;
    using verlet_list_type = mjolnir::VerletList<traits_type, potential_type, false>;

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
    std::iota(participants.begin(), participants.end(), 0);

    dummy_potential<real_type> pot(cutoff, participants);

    mjolnir::System<traits_type> sys(N, boundary_type(coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    sys.topology().construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(
            mjolnir::make_unique<verlet_list_type>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : participants)
    {
        const auto partners = vlist.partners(i);
        for(std::size_t j=0; j<N; ++j)
        {
            if(std::find_if(partners.begin(), partners.end(),
                [=](const neighbor_type& elem) -> bool {return elem.index == j;}
                        ) == partners.end())
            {
                if(i == j){continue;}

                // should be enough distant (>= threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist >= threshold);
            }
            else
            {
                BOOST_TEST(i != j);

                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist < threshold);
            }
        }
    }

    // check parameter_type.
    for(const auto i : participants)
    {
        for(const auto& p_j : vlist.partners(i))
        {
            const std::size_t j = p_j.index;
            BOOST_TEST(p_j.parameter().first  == pot.prepare_params(i, j).first);
            BOOST_TEST(p_j.parameter().second == pot.prepare_params(i, j).second);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_VerletList_PeriodicBoundary_partial_interchange)
{
    mjolnir::LoggerManager::set_default_logger("test_verlet_list.log");
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = dummy_potential<real_type>;
    using verlet_list_type = mjolnir::VerletList<traits_type, potential_type, false>;

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

    dummy_potential<real_type> pot(cutoff, participants);

    mjolnir::System<traits_type> sys(N, boundary_type(coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    sys.topology().construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(
        mjolnir::make_unique<verlet_list_type>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : participants)
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

        for(std::size_t j=0; j<N; ++j)
        {
            if(std::find_if(partners.begin(), partners.end(),
                [=](const neighbor_type& elem) -> bool {return elem.index == j;}
                        ) == partners.end())
            {
                if(i == j){continue;}

                // should be enough distant (>= threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                const bool enough_distant = dist >= threshold;

                // or not a participant
                const bool is_participant = pot.participants().end() != std::find(
                    pot.participants().begin(), pot.participants().end(), j);

                const bool is_ok = enough_distant || (!is_participant);
                BOOST_TEST(is_ok);
            }
            else
            {
                BOOST_TEST(i != j);

                // should be a participant
                const bool found = pot.participants().end() != std::find(
                    pot.participants().begin(), pot.participants().end(), j);
                BOOST_TEST(found);

                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_VerletList_PeriodicBoundary_partial_2_interchange)
{
    mjolnir::LoggerManager::set_default_logger("test_verlet_list.log");
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>;
    using real_type       = typename traits_type::real_type;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using potential_type  = dummy_potential<real_type>;
    using verlet_list_type = mjolnir::VerletList<traits_type, potential_type, false>;

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

    dummy_potential<real_type> pot(cutoff, participants);

    mjolnir::System<traits_type> sys(N, boundary_type(coordinate_type(0.0, 0.0, 0.0), coordinate_type(L, L, L)));

    std::mt19937 mt(123456789);
    for(std::size_t i=0; i < N; ++i)
    {
        sys.at(i).mass     = 1.0;
        sys.at(i).position = distribute_particle(mt, L);
    }
    sys.topology().construct_molecules();

    mjolnir::SpatialPartition<traits_type, potential_type> vlist(
        mjolnir::make_unique<verlet_list_type>(margin));

    using neighbor_type = typename decltype(vlist)::neighbor_type;

    BOOST_TEST(!vlist.valid());

    vlist.initialize(sys, pot);
    vlist.make(sys, pot);
    BOOST_TEST(vlist.valid());

    for(const auto i : participants)
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

        for(std::size_t j=0; j<N; ++j)
        {
            if(std::find_if(partners.begin(), partners.end(),
                [=](const neighbor_type& elem) -> bool {return elem.index == j;}
                        ) == partners.end())
            {
                if(i == j){continue;}

                // should be enough distant (>= threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                const bool enough_distant = dist >= threshold;

                // or not a participant
                const bool is_participant = pot.participants().end() != std::find(
                    pot.participants().begin(), pot.participants().end(), j);

                const bool is_ok = enough_distant || (!is_participant);
                BOOST_TEST(is_ok);
            }
            else
            {
                BOOST_TEST(i != j);

                // should be a participant
                const bool found = pot.participants().end() != std::find(
                    pot.participants().begin(), pot.participants().end(), j);
                BOOST_TEST(found);

                // should be enough close (< threshold)
                const auto dist = mjolnir::math::length(sys.adjust_direction(
                            sys.position(j) - sys.position(i)));
                BOOST_TEST(dist < threshold);
            }
        }
    }
}
