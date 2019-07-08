#define BOOST_TEST_MODULE "test_exclusion_list"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/empty.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/potential/global/IgnoreMolecule.hpp>
#include <limits>

#include <mjolnir/core/ExclusionList.hpp>
#include <random>
#include <cmath>

struct dummy_potential
{
  public:
    using real_type = double;
    using parameter_type = mjolnir::empty_t;

    using topology_type        = mjolnir::Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = mjolnir::IgnoreMolecule<molecule_id_type>;

  public:

    dummy_potential(const std::map<connection_kind_type, std::size_t>& exclusions,
                    ignore_molecule_type ignore_mol)
        : ignore_molecule_(std::move(ignore_mol)),
          ignore_within_  (exclusions.begin(), exclusions.end())
    {}

    parameter_type prepare_params(std::size_t, std::size_t) const noexcept {return {};}

    // forwarding functions for clarity...
    real_type potential(const std::size_t, const std::size_t, const real_type) const noexcept
    {
        return 0.0;
    }
    real_type derivative(const std::size_t, const std::size_t, const real_type) const noexcept
    {
        return 0.0;
    }

    real_type potential(const real_type, const parameter_type&) const noexcept
    {
        return 0.0;
    }
    real_type derivative(const real_type, const parameter_type&) const noexcept
    {
        return 0.0;
    }

    template<typename traitsT>
    void initialize(const mjolnir::System<traitsT>&) noexcept {return;}

    // nothing to be done if system parameter (e.g. temperature) changes
    template<typename traitsT>
    void update(const mjolnir::System<traitsT>&) const noexcept {return;}

    real_type max_cutoff_length() const
    {
        return std::numeric_limits<real_type>::max();
    }

    std::vector<std::pair<connection_kind_type, std::size_t>>
    ignore_within() const {return ignore_within_;}

    bool is_ignored_molecule(
            const molecule_id_type& i, const molecule_id_type& j) const noexcept
    {
        return ignore_molecule_.is_ignored(i, j);
    }

    const char* name() const noexcept {return "dummy";}

  private:
    ignore_molecule_type ignore_molecule_;
    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within_;
};

BOOST_AUTO_TEST_CASE(ExclusionList_noignore)
{
    mjolnir::LoggerManager::set_default_logger("test_ExclusionList");
    using traits_type          = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using boundary_type        = traits_type::boundary_type;
    using topology_type        = mjolnir::Topology;
    using molecule_id_type     = topology_type::molecule_id_type;
    using ignore_molecule_type = mjolnir::IgnoreMolecule<molecule_id_type>;

    const dummy_potential pot({}, ignore_molecule_type("Nothing"));
    constexpr std::size_t N = 10;

    // no topology
    {
        mjolnir::System<traits_type> sys(N, boundary_type{});
        for(std::size_t i=0; i<N; ++i)
        {
            sys.name(i)  = "X";
            sys.group(i) = "none";
        }
        sys.topology().construct_molecules();

        mjolnir::ExclusionList exl;
        exl.make(sys, pot);

        for(std::size_t i=0; i<10; ++i)
        {
            for(std::size_t j=0; j<10; ++j)
            {
                if(i == j)
                {
                    BOOST_TEST(exl.is_excluded(i, j));
                }
                else
                {
                    BOOST_TEST(!exl.is_excluded(i, j));
                }
            }
        }
    }

    // add topology (no effect)

    {
        mjolnir::System<traits_type> sys(N, boundary_type{});
        for(std::size_t i=0; i<N; ++i)
        {
            sys.name(i)  = "X";
            sys.group(i) = "none";
        }

        for(std::size_t i=1; i<N; ++i)
        {
            sys.topology().add_connection(i-1, i, "bond");
        }
        sys.topology().add_connection(0, N-1, "contact");

        sys.topology().construct_molecules();

        mjolnir::ExclusionList exl;
        exl.make(sys, pot);

        for(std::size_t i=0; i<10; ++i)
        {
            for(std::size_t j=0; j<10; ++j)
            {
                if(i == j)
                {
                    BOOST_TEST(exl.is_excluded(i, j));
                }
                else
                {
                    BOOST_TEST(!exl.is_excluded(i, j));
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(ExclusionList_topology_dependent)
{
    mjolnir::LoggerManager::set_default_logger("test_ExclusionList");
    using traits_type          = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using boundary_type        = traits_type::boundary_type;
    using topology_type        = mjolnir::Topology;
    using molecule_id_type     = topology_type::molecule_id_type;
    using ignore_molecule_type = mjolnir::IgnoreMolecule<molecule_id_type>;

    const dummy_potential pot({{"bond", 3}, {"contact", 1}},
            ignore_molecule_type("Nothing"));

    // no topology
    {
        mjolnir::System<traits_type> sys(10, boundary_type{});
        for(std::size_t i=0; i<10; ++i)
        {
            sys.name(i)  = "X";
            sys.group(i) = "none";
        }
        sys.topology().construct_molecules();

        mjolnir::ExclusionList exl;
        exl.make(sys, pot);

        for(std::size_t i=0; i<10; ++i)
        {
            for(std::size_t j=0; j<10; ++j)
            {
                if(i == j)
                {
                    BOOST_TEST(exl.is_excluded(i, j));
                }
                else
                {
                    BOOST_TEST(!exl.is_excluded(i, j));
                }
            }
        }
    }

    // add topology
    {
        mjolnir::System<traits_type> sys(10, boundary_type{});
        for(std::size_t i=0; i<10; ++i)
        {
            sys.name(i)  = "X";
            sys.group(i) = "none";
        }

        for(std::size_t i=1; i<10; ++i)
        {
            sys.topology().add_connection(i-1, i, "bond");
        }
        sys.topology().add_connection(0, 9, "contact");

        sys.topology().construct_molecules();

        mjolnir::ExclusionList exl;
        exl.make(sys, pot);

        for(std::int32_t i=0; i<10; ++i)
        {
            for(std::int32_t j=0; j<10; ++j)
            {
                if(std::abs(i - j) <= 3) // bond * 3
                {
                    BOOST_TEST(exl.is_excluded(i, j));
                }
                else if((i==0 && j==9) || (i==9 && j == 0)) // contact
                {
                    BOOST_TEST(exl.is_excluded(i, j));
                }
                else
                {
                    BOOST_TEST(!exl.is_excluded(i, j));
                }
            }
        }
    }
    // add branch
    {
        //        0 -- 1  |
        //       /        |
        // 3 -- 2         |
        //       \        |
        //        4 -- 5  |
        //       /        |
        // 7 -- 6         |
        //       \        |
        //        8 -- 9  |

        mjolnir::System<traits_type> sys(10, boundary_type{});
        for(std::size_t i=0; i<10; ++i)
        {
            sys.name(i)  = "X";
            sys.group(i) = "none";
        }

        sys.topology().add_connection(0, 1, "bond");
        sys.topology().add_connection(0, 2, "bond");
        sys.topology().add_connection(2, 3, "bond");
        sys.topology().add_connection(2, 4, "bond");
        sys.topology().add_connection(4, 5, "bond");
        sys.topology().add_connection(4, 6, "bond");
        sys.topology().add_connection(6, 7, "bond");
        sys.topology().add_connection(6, 8, "bond");
        sys.topology().add_connection(8, 9, "bond");

        sys.topology().construct_molecules();

        mjolnir::ExclusionList exl;
        exl.make(sys, pot);

        BOOST_TEST( exl.is_excluded(0, 0));
        BOOST_TEST( exl.is_excluded(0, 1));
        BOOST_TEST( exl.is_excluded(0, 2));
        BOOST_TEST( exl.is_excluded(0, 3));
        BOOST_TEST( exl.is_excluded(0, 4));
        BOOST_TEST( exl.is_excluded(0, 5));
        BOOST_TEST( exl.is_excluded(0, 6));
        BOOST_TEST(!exl.is_excluded(0, 7));
        BOOST_TEST(!exl.is_excluded(0, 8));
        BOOST_TEST(!exl.is_excluded(0, 9));

        BOOST_TEST( exl.is_excluded(1, 0));
        BOOST_TEST( exl.is_excluded(1, 1));
        BOOST_TEST( exl.is_excluded(1, 2));
        BOOST_TEST( exl.is_excluded(1, 3));
        BOOST_TEST( exl.is_excluded(1, 4));
        BOOST_TEST(!exl.is_excluded(1, 5));
        BOOST_TEST(!exl.is_excluded(1, 6));
        BOOST_TEST(!exl.is_excluded(1, 7));
        BOOST_TEST(!exl.is_excluded(1, 8));
        BOOST_TEST(!exl.is_excluded(1, 9));

        BOOST_TEST( exl.is_excluded(2, 0));
        BOOST_TEST( exl.is_excluded(2, 1));
        BOOST_TEST( exl.is_excluded(2, 2));
        BOOST_TEST( exl.is_excluded(2, 3));
        BOOST_TEST( exl.is_excluded(2, 4));
        BOOST_TEST( exl.is_excluded(2, 5));
        BOOST_TEST( exl.is_excluded(2, 6));
        BOOST_TEST( exl.is_excluded(2, 7));
        BOOST_TEST( exl.is_excluded(2, 8));
        BOOST_TEST(!exl.is_excluded(2, 9));

        BOOST_TEST( exl.is_excluded(3, 0));
        BOOST_TEST( exl.is_excluded(3, 1));
        BOOST_TEST( exl.is_excluded(3, 2));
        BOOST_TEST( exl.is_excluded(3, 3));
        BOOST_TEST( exl.is_excluded(3, 4));
        BOOST_TEST( exl.is_excluded(3, 5));
        BOOST_TEST( exl.is_excluded(3, 6));
        BOOST_TEST(!exl.is_excluded(3, 7));
        BOOST_TEST(!exl.is_excluded(3, 8));
        BOOST_TEST(!exl.is_excluded(3, 9));

        BOOST_TEST( exl.is_excluded(4, 0));
        BOOST_TEST( exl.is_excluded(4, 1));
        BOOST_TEST( exl.is_excluded(4, 2));
        BOOST_TEST( exl.is_excluded(4, 3));
        BOOST_TEST( exl.is_excluded(4, 4));
        BOOST_TEST( exl.is_excluded(4, 5));
        BOOST_TEST( exl.is_excluded(4, 6));
        BOOST_TEST( exl.is_excluded(4, 7));
        BOOST_TEST( exl.is_excluded(4, 8));
        BOOST_TEST( exl.is_excluded(4, 9));

        //        0 -- 1  |
        //       /        |
        // 3 -- 2         |
        //       \        |
        //        4 -- 5  |
        //       /        |
        // 7 -- 6         |
        //       \        |
        //        8 -- 9  |

        BOOST_TEST( exl.is_excluded(5, 0));
        BOOST_TEST(!exl.is_excluded(5, 1));
        BOOST_TEST( exl.is_excluded(5, 2));
        BOOST_TEST( exl.is_excluded(5, 3));
        BOOST_TEST( exl.is_excluded(5, 4));
        BOOST_TEST( exl.is_excluded(5, 5));
        BOOST_TEST( exl.is_excluded(5, 6));
        BOOST_TEST( exl.is_excluded(5, 7));
        BOOST_TEST( exl.is_excluded(5, 8));
        BOOST_TEST(!exl.is_excluded(5, 9));

        BOOST_TEST( exl.is_excluded(6, 0));
        BOOST_TEST(!exl.is_excluded(6, 1));
        BOOST_TEST( exl.is_excluded(6, 2));
        BOOST_TEST( exl.is_excluded(6, 3));
        BOOST_TEST( exl.is_excluded(6, 4));
        BOOST_TEST( exl.is_excluded(6, 5));
        BOOST_TEST( exl.is_excluded(6, 6));
        BOOST_TEST( exl.is_excluded(6, 7));
        BOOST_TEST( exl.is_excluded(6, 8));
        BOOST_TEST( exl.is_excluded(6, 9));

        BOOST_TEST(!exl.is_excluded(7, 0));
        BOOST_TEST(!exl.is_excluded(7, 1));
        BOOST_TEST( exl.is_excluded(7, 2));
        BOOST_TEST(!exl.is_excluded(7, 3));
        BOOST_TEST( exl.is_excluded(7, 4));
        BOOST_TEST( exl.is_excluded(7, 5));
        BOOST_TEST( exl.is_excluded(7, 6));
        BOOST_TEST( exl.is_excluded(7, 7));
        BOOST_TEST( exl.is_excluded(7, 8));
        BOOST_TEST( exl.is_excluded(7, 9));

        BOOST_TEST(!exl.is_excluded(8, 0));
        BOOST_TEST(!exl.is_excluded(8, 1));
        BOOST_TEST( exl.is_excluded(8, 2));
        BOOST_TEST(!exl.is_excluded(8, 3));
        BOOST_TEST( exl.is_excluded(8, 4));
        BOOST_TEST( exl.is_excluded(8, 5));
        BOOST_TEST( exl.is_excluded(8, 6));
        BOOST_TEST( exl.is_excluded(8, 7));
        BOOST_TEST( exl.is_excluded(8, 8));
        BOOST_TEST( exl.is_excluded(8, 9));

        BOOST_TEST(!exl.is_excluded(9, 0));
        BOOST_TEST(!exl.is_excluded(9, 1));
        BOOST_TEST(!exl.is_excluded(9, 2));
        BOOST_TEST(!exl.is_excluded(9, 3));
        BOOST_TEST( exl.is_excluded(9, 4));
        BOOST_TEST(!exl.is_excluded(9, 5));
        BOOST_TEST( exl.is_excluded(9, 6));
        BOOST_TEST( exl.is_excluded(9, 7));
        BOOST_TEST( exl.is_excluded(9, 8));
        BOOST_TEST( exl.is_excluded(9, 9));
    }
}

BOOST_AUTO_TEST_CASE(ExclusionList_molecule_dependent)
{
    mjolnir::LoggerManager::set_default_logger("test_ExclusionList");
    using traits_type          = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using boundary_type        = traits_type::boundary_type;
    using topology_type        = mjolnir::Topology;
    using molecule_id_type     = topology_type::molecule_id_type;
    using ignore_molecule_type = mjolnir::IgnoreMolecule<molecule_id_type>;

    {
        // there are no interaction between particles in the same molecules
        const dummy_potential pot({{"bond", 1}},
                ignore_molecule_type("Self"));

        {
            mjolnir::System<traits_type> sys(10, boundary_type{});
            for(std::size_t i=0; i<10; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "none";
            }

            for(std::size_t i=1; i<10; ++i)
            {
                sys.topology().add_connection(i-1, i, "bond");
            }
            sys.topology().construct_molecules();

            mjolnir::ExclusionList exl;
            exl.make(sys, pot);

            for(std::int32_t i=0; i<10; ++i)
            {
                for(std::int32_t j=0; j<10; ++j)
                {
                    BOOST_TEST(exl.is_excluded(i, j));
                }
            }
        }
        {
            mjolnir::System<traits_type> sys(10, boundary_type{});
            for(std::size_t i=0; i<10; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "none";
            }

            for(std::size_t i=1; i<10; ++i)
            {
                sys.topology().add_connection(i-1, i, "bond");
            }
            sys.topology().erase_connection(4, 5, "bond");

            sys.topology().construct_molecules();

            mjolnir::ExclusionList exl;
            exl.make(sys, pot);

            for(std::int32_t i=0; i<10; ++i)
            {
                for(std::int32_t j=0; j<10; ++j)
                {
                    if((i <= 4 && 5 <= j) || (j <= 4 && 5 <= i))
                    {
                        BOOST_TEST(!exl.is_excluded(i, j));
                    }
                    else
                    {
                        BOOST_TEST( exl.is_excluded(i, j));
                    }
                }
            }
        }
    }

    {
        // there are no interaction between particles in different molecules
        const dummy_potential pot({{"bond", 1}},
                ignore_molecule_type("Others"));

        {
            mjolnir::System<traits_type> sys(10, boundary_type{});
            for(std::size_t i=0; i<10; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "none";
            }

            for(std::size_t i=1; i<10; ++i)
            {
                sys.topology().add_connection(i-1, i, "bond");
            }
            sys.topology().construct_molecules();

            mjolnir::ExclusionList exl;
            exl.make(sys, pot);

            for(std::int32_t i=0; i<10; ++i)
            {
                for(std::int32_t j=0; j<10; ++j)
                {
                    if(std::abs(i - j) <= 1)
                    {
                        BOOST_TEST( exl.is_excluded(i, j));
                    }
                    else
                    {
                        BOOST_TEST(!exl.is_excluded(i, j));
                    }
                }
            }
        }
        {
            mjolnir::System<traits_type> sys(10, boundary_type{});
            for(std::size_t i=0; i<10; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "none";
            }

            for(std::size_t i=1; i<10; ++i)
            {
                sys.topology().add_connection(i-1, i, "bond");
            }
            sys.topology().erase_connection(4, 5, "bond");

            sys.topology().construct_molecules();

            BOOST_TEST(sys.topology().number_of_molecules() == 2);

            mjolnir::ExclusionList exl;
            exl.make(sys, pot);

            for(std::int32_t i=0; i<10; ++i)
            {
                for(std::int32_t j=0; j<10; ++j)
                {
                    if(sys.topology().molecule_of(i) != sys.topology().molecule_of(j))
                    {
                        BOOST_TEST( exl.is_excluded(i, j));
                    }
                    else
                    {
                        if(std::abs(i - j) <= 1)
                        {
                            BOOST_TEST( exl.is_excluded(i, j));
                        }
                        else
                        {
                            BOOST_TEST(!exl.is_excluded(i, j));
                        }
                    }
                }
            }
        }
    }
}
