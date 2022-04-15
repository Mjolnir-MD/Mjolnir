#define BOOST_TEST_MODULE "test_read_uniform_cubic_pan_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_global_potential.hpp>

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_cubic_pan_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_cubic_pan.log");
    using real_type = T;

    // a dummy system for testing `initialize` method
    using traits_type   = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::UniformCubicPanPotential<real_type>;
    using boundary_type = typename traits_type::boundary_type;
    mjolnir::System<traits_type> sys(10, boundary_type{});
    mjolnir::Topology          topol(10);
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformCubicPan"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            epsilon = 1.5
            v0      = 5.0
            range   = 5.0
        )"_toml;

        auto pot_para = mjolnir::read_uniform_cubic_pan_potential<traits_type>(v);
        const auto& pot = pot_para.first;
        auto& para = dynamic_cast<mjolnir::EmptyCombinationRule<traits_type, potential_type>&>(pot_para.second.ref());

        BOOST_TEST(para.participants().empty());
        para.initialize(sys, topol, pot);

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(pot.epsilon()           == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(pot.v0()                == real_type(5.0), tolerance<real_type>());
        BOOST_TEST(pot.interaction_range() == real_type(5.0), tolerance<real_type>());

        BOOST_TEST(para.participants().size() == 10u);
        BOOST_TEST(para.participants().at(0) == 0u);
        BOOST_TEST(para.participants().at(1) == 1u);
        BOOST_TEST(para.participants().at(2) == 2u);
        BOOST_TEST(para.participants().at(3) == 3u);
        BOOST_TEST(para.participants().at(4) == 4u);
        BOOST_TEST(para.participants().at(5) == 5u);
        BOOST_TEST(para.participants().at(6) == 6u);
        BOOST_TEST(para.participants().at(7) == 7u);
        BOOST_TEST(para.participants().at(8) == 8u);
        BOOST_TEST(para.participants().at(9) == 9u);
    }
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformCubicPan"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            "ε"   = 1.5
            v0    = 5.0
            range = 5.0
        )"_toml;

        auto pot_para = mjolnir::read_uniform_cubic_pan_potential<traits_type>(v);
        const auto& pot = pot_para.first;
        auto& para = dynamic_cast<mjolnir::EmptyCombinationRule<traits_type, potential_type>&>(pot_para.second.ref());

        BOOST_TEST(para.participants().empty());
        para.initialize(sys, topol, pot);

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(pot.epsilon()           == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(pot.v0()                == real_type(5.0), tolerance<real_type>());
        BOOST_TEST(pot.interaction_range() == real_type(5.0), tolerance<real_type>());

        BOOST_TEST(para.participants().size() == 10u);
        BOOST_TEST(para.participants().at(0) == 0u);
        BOOST_TEST(para.participants().at(1) == 1u);
        BOOST_TEST(para.participants().at(2) == 2u);
        BOOST_TEST(para.participants().at(3) == 3u);
        BOOST_TEST(para.participants().at(4) == 4u);
        BOOST_TEST(para.participants().at(5) == 5u);
        BOOST_TEST(para.participants().at(6) == 6u);
        BOOST_TEST(para.participants().at(7) == 7u);
        BOOST_TEST(para.participants().at(8) == 8u);
        BOOST_TEST(para.participants().at(9) == 9u);
    }
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformCubicPan"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            "ε"   = 1.5
            v0    = 5.0
            range = 5.0
            parameters = [
                {index = 1},
                {index = 2},
                {index = 3},
                {index = 4},
                {index = 5},
            ]
        )"_toml;

        auto pot_para = mjolnir::read_uniform_cubic_pan_potential<traits_type>(v);
        const auto& pot = pot_para.first;
        auto& para = dynamic_cast<mjolnir::EmptyCombinationRule<traits_type, potential_type>&>(pot_para.second.ref());

        para.initialize(sys, topol, pot);

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(pot.epsilon()           == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(pot.v0()                == real_type(5.0), tolerance<real_type>());
        BOOST_TEST(pot.interaction_range() == real_type(5.0), tolerance<real_type>());

        BOOST_TEST(para.participants().size() == 5u);
        BOOST_TEST(para.participants().at(0)  == 1u);
        BOOST_TEST(para.participants().at(1)  == 2u);
        BOOST_TEST(para.participants().at(2)  == 3u);
        BOOST_TEST(para.participants().at(3)  == 4u);
        BOOST_TEST(para.participants().at(4)  == 5u);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_cubic_pan_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_cubic_pan.log");
    using real_type = T;

    // a dummy system for testing `initialize` method
    using traits_type   = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::UniformCubicPanPotential<real_type>;
    using boundary_type = typename traits_type::boundary_type;
    mjolnir::System<traits_type> sys(10, boundary_type{});
    mjolnir::Topology          topol(10);

    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformCubicPan"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            "ε"   = 1.5
            v0    = 5.0
            range = 5.0
            env.three = 3
            parameters = [
                {index = 1},
                {index = 2},
                {index = "three"},
                {index = 4},
                {index = 5},
            ]
        )"_toml;

        auto pot_para = mjolnir::read_uniform_cubic_pan_potential<traits_type>(v);
        const auto& pot = pot_para.first;
        auto& para = dynamic_cast<mjolnir::EmptyCombinationRule<traits_type, potential_type>&>(pot_para.second.ref());

        para.initialize(sys, topol, pot);

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);
        BOOST_TEST(pot.epsilon()           == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(pot.v0()                == real_type(5.0), tolerance<real_type>());
        BOOST_TEST(pot.interaction_range() == real_type(5.0), tolerance<real_type>());

        BOOST_TEST(para.participants().size() == 5u);
        BOOST_TEST(para.participants().at(0)  == 1u);
        BOOST_TEST(para.participants().at(1)  == 2u);
        BOOST_TEST(para.participants().at(2)  == 3u);
        BOOST_TEST(para.participants().at(3)  == 4u);
        BOOST_TEST(para.participants().at(4)  == 5u);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_cubic_pan_ignore_self, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_cubic_pan.log");
    using real_type = T;

    // a dummy system for testing `initialize` method
    using traits_type   = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::UniformCubicPanPotential<real_type>;
    using boundary_type = typename traits_type::boundary_type;
    mjolnir::System<traits_type> sys(10, boundary_type{});
    mjolnir::Topology            topol(10);

    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformCubicPan"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Self"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            epsilon = 1.5
            v0      = 5.0
            range   = 5.0
            parameters = []
        )"_toml;

        auto pot_para = mjolnir::read_uniform_cubic_pan_potential<traits_type>(v);
        const auto& pot = pot_para.first;
        auto& para = dynamic_cast<mjolnir::EmptyCombinationRule<traits_type, potential_type>&>(pot_para.second.ref());

        BOOST_TEST(para.participants().empty());
        para.initialize(sys, topol, pot);

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST( para.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST( para.exclusion_list().is_ignored_molecule(1, 1));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_cubic_pan_ignore_others, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_cubic_pan.log");
    using real_type = T;

    // a dummy system for testing `initialize` method
    using traits_type   = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::UniformCubicPanPotential<real_type>;
    using boundary_type = typename traits_type::boundary_type;
    mjolnir::System<traits_type> sys(10, boundary_type{});
    mjolnir::Topology            topol(10);

    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformCubicPan"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Others"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            epsilon = 1.5
            v0      = 5.0
            range   = 5.0
            parameters = []
        )"_toml;

        auto pot_para = mjolnir::read_uniform_cubic_pan_potential<traits_type>(v);
        const auto& pot = pot_para.first;
        auto& para = dynamic_cast<mjolnir::EmptyCombinationRule<traits_type, potential_type>&>(pot_para.second.ref());

        BOOST_TEST(para.participants().empty());
        para.initialize(sys, topol, pot);

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST( para.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(1, 1));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_uniform_cubic_pan_ignore_group, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_uniform_cubic_pan.log");
    using real_type = T;

    // a dummy system for testing `initialize` method
    using traits_type   = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::UniformCubicPanPotential<real_type>;
    using boundary_type = typename traits_type::boundary_type;
    mjolnir::System<traits_type> sys(10, boundary_type{});
    mjolnir::Topology            topol(10);

    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformCubicPan"
            spatial_partition.type  = "CellList"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            ignore.group.inter      = [
                ["protein1", "protein2"], # between these
                ["protein1", "protein3"],
            ]
            epsilon = 1.5
            v0      = 5.0
            range   = 5.0
            parameters = []
        )"_toml;

        auto pot_para = mjolnir::read_uniform_cubic_pan_potential<traits_type>(v);
        const auto& pot = pot_para.first;
        auto& para = dynamic_cast<mjolnir::EmptyCombinationRule<traits_type, potential_type>&>(pot_para.second.ref());

        BOOST_TEST(para.participants().empty());
        para.initialize(sys, topol, pot);

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(1, 1));

        BOOST_TEST(!para.exclusion_list().is_ignored_group("protein1", "protein1"));
        BOOST_TEST( para.exclusion_list().is_ignored_group("protein1", "protein2"));
        BOOST_TEST( para.exclusion_list().is_ignored_group("protein1", "protein3"));

        BOOST_TEST( para.exclusion_list().is_ignored_group("protein2", "protein1"));
        BOOST_TEST(!para.exclusion_list().is_ignored_group("protein2", "protein2"));
        BOOST_TEST(!para.exclusion_list().is_ignored_group("protein2", "protein3"));

        BOOST_TEST( para.exclusion_list().is_ignored_group("protein3", "protein1"));
        BOOST_TEST(!para.exclusion_list().is_ignored_group("protein3", "protein2"));
        BOOST_TEST(!para.exclusion_list().is_ignored_group("protein3", "protein3"));
    }
}
