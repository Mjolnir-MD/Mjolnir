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
#include <mjolnir/potential/global/IgnoreGroup.hpp>
#include <mjolnir/potential/global/IgnoreMolecule.hpp>
#include <limits>

#include <mjolnir/core/ExclusionList.hpp>
#include <random>
#include <cmath>

BOOST_AUTO_TEST_CASE(ExclusionList_noignore)
{
    mjolnir::LoggerManager::set_default_logger("test_ExclusionList");
    using traits_type          = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using boundary_type        = traits_type::boundary_type;
    using topology_type        = mjolnir::Topology;
    using molecule_id_type     = topology_type::molecule_id_type;
    using group_id_type        = topology_type::group_id_type;
    using ignore_molecule_type = mjolnir::IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = mjolnir::IgnoreGroup   <group_id_type>;

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

        mjolnir::ExclusionList exl({},
                ignore_molecule_type("Nothing"), ignore_group_type({}));
        exl.make(sys);

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

        mjolnir::ExclusionList exl({},
                ignore_molecule_type("Nothing"), ignore_group_type({}));
        exl.make(sys);

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
    using group_id_type        = topology_type::group_id_type;
    using ignore_molecule_type = mjolnir::IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = mjolnir::IgnoreGroup   <group_id_type>;

    // no topology
    {
        mjolnir::System<traits_type> sys(10, boundary_type{});
        for(std::size_t i=0; i<10; ++i)
        {
            sys.name(i)  = "X";
            sys.group(i) = "none";
        }
        sys.topology().construct_molecules();

        mjolnir::ExclusionList exl({{"bond", 3}, {"contact", 1}},
                                   ignore_molecule_type("Nothing"),
                                   ignore_group_type({}));
        exl.make(sys);

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

        mjolnir::ExclusionList exl({{"bond", 3}, {"contact", 1}},
                                   ignore_molecule_type("Nothing"),
                                   ignore_group_type({}));
        exl.make(sys);

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

        mjolnir::ExclusionList exl({{"bond", 3}, {"contact", 1}},
                                   ignore_molecule_type("Nothing"),
                                   ignore_group_type({}));
        exl.make(sys);

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
    using group_id_type        = topology_type::group_id_type;
    using ignore_molecule_type = mjolnir::IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = mjolnir::IgnoreGroup   <group_id_type>;

    {
        // there are no interaction between particles in the same molecules
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

            mjolnir::ExclusionList exl({{"bond", 1}},
                ignore_molecule_type("Self"),
                ignore_group_type({}));
            exl.make(sys);

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

            mjolnir::ExclusionList exl({{"bond", 1}},
                ignore_molecule_type("Self"),
                ignore_group_type({}));
            exl.make(sys);

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

            mjolnir::ExclusionList exl({{"bond", 1}},
                ignore_molecule_type("Others"),
                ignore_group_type({}));
            exl.make(sys);

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

            BOOST_TEST(sys.topology().number_of_molecules() == 2u);

            mjolnir::ExclusionList exl({{"bond", 1}},
                ignore_molecule_type("Others"),
                ignore_group_type({}));
            exl.make(sys);

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

BOOST_AUTO_TEST_CASE(ExclusionList_gruop_dependent)
{
    mjolnir::LoggerManager::set_default_logger("test_ExclusionList");
    using traits_type          = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using boundary_type        = traits_type::boundary_type;
    using topology_type        = mjolnir::Topology;
    using molecule_id_type     = topology_type::molecule_id_type;
    using group_id_type        = topology_type::group_id_type;
    using ignore_molecule_type = mjolnir::IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = mjolnir::IgnoreGroup   <group_id_type>;

    {
        {
            mjolnir::System<traits_type> sys(10, boundary_type{});
            for(std::size_t i=0; i<5; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "protein1";
            }
            for(std::size_t i=5; i<10; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "protein2";
            }

            for(std::size_t i=1; i<10; ++i)
            {
                sys.topology().add_connection(i-1, i, "bond");
            }
            sys.topology().erase_connection(4, 5, "bond");
            sys.topology().construct_molecules();

            mjolnir::ExclusionList exl({{"bond", 1}},
                ignore_molecule_type("Nothing"),
                ignore_group_type({
                    // interact only if groups are different (intra-groups are ignored)
                    {"protein1", {"protein1"}},
                    {"protein2", {"protein2"}}
                }));
            exl.make(sys);

            for(std::int32_t i=0; i<10; ++i)
            {
                for(std::int32_t j=0; j<10; ++j)
                {
                    if(sys.group(i) != sys.group(j))
                    {
                        BOOST_TEST(!exl.is_excluded(i, j));
                    }
                    else
                    {
                        BOOST_TEST(exl.is_excluded(i, j));
                    }
                }
            }
        }
    }

    {
        {
            mjolnir::System<traits_type> sys(10, boundary_type{});
            for(std::size_t i=0; i<5; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "protein1";
            }
            for(std::size_t i=5; i<10; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "protein2";
            }

            for(std::size_t i=1; i<10; ++i)
            {
                sys.topology().add_connection(i-1, i, "bond");
            }
            sys.topology().erase_connection(4, 5, "bond");
            sys.topology().construct_molecules();

            mjolnir::ExclusionList exl({{"bond", 1}},
                ignore_molecule_type("Nothing"),
                ignore_group_type({
                    // interact only if groups are different (inter-groups are ignored)
                    {"protein1", {"protein2"}},
                    {"protein2", {"protein1"}}
                }));
            exl.make(sys);

            for(std::int32_t i=0; i<10; ++i)
            {
                for(std::int32_t j=0; j<10; ++j)
                {
                    if(sys.group(i) != sys.group(j))
                    {
                        BOOST_TEST(exl.is_excluded(i, j));
                    }
                    else // the same group
                    {
                        if(sys.topology().has_connection(i, j, "bond"))
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
    }

    {
        {
            mjolnir::System<traits_type> sys(10, boundary_type{});
            for(std::size_t i=0; i<5; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "protein1";
            }
            for(std::size_t i=5; i<10; ++i)
            {
                sys.name(i)  = "X";
                sys.group(i) = "protein2";
            }

            for(std::size_t i=1; i<10; ++i)
            {
                sys.topology().add_connection(i-1, i, "bond");
            }
            sys.topology().erase_connection(4, 5, "bond");
            sys.topology().construct_molecules();

            mjolnir::ExclusionList exl({{"bond", 1}},
                ignore_molecule_type("Nothing"),
                ignore_group_type({
                    // everything will be ignored
                    {"protein1", {"protein1", "protein2"}},
                    {"protein2", {"protein1", "protein2"}}
                }));
            exl.make(sys);

            for(std::int32_t i=0; i<10; ++i)
            {
                for(std::int32_t j=0; j<10; ++j)
                {
                    BOOST_TEST(exl.is_excluded(i, j));
                }
            }
        }
    }
}
