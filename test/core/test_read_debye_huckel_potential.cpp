#define BOOST_TEST_MODULE "test_read_debye_huckel_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_debye_huckel_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_debye_huckel.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction                     = "Pair"
            potential                       = "DebyeHuckel"
            spatial_partition.type          = "Nothing"
            ignore.molecule                 = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = [
                {index =   0, charge =   1.0},
                {index =   1, charge =  -1.0},
                {index =   3, charge =   0.3},
                {index =   5, charge =   0.5},
                {index =   7, charge =   0.7},
                {index = 100, charge = 100.0},
            ]
        )"_toml;

        const auto pot = mjolnir::read_debye_huckel_potential<traits_type>(v);

        const auto ignore_within = pot.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(pot.participants().size() ==   6u);
        BOOST_TEST(pot.participants().at(0)  ==   0u);
        BOOST_TEST(pot.participants().at(1)  ==   1u);
        BOOST_TEST(pot.participants().at(2)  ==   3u);
        BOOST_TEST(pot.participants().at(3)  ==   5u);
        BOOST_TEST(pot.participants().at(4)  ==   7u);
        BOOST_TEST(pot.participants().at(5)  == 100u);

        BOOST_TEST(pot.charges().at(  0)  == real_type(  1.0), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(  1)  == real_type( -1.0), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(  3)  == real_type(  0.3), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(  5)  == real_type(  0.5), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(  7)  == real_type(  0.7), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(100)  == real_type(100.0), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_debye_huckel_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_debye_huckel.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction                     = "Pair"
            potential                       = "DebyeHuckel"
            spatial_partition.type          = "Nothing"
            ignore.molecule                 = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            env.positive = 1.0
            env.negative = -1.0
            env.toomuch  = 100.0
            parameters = [
                {index =   0, charge = "positive"},
                {index =   1, charge = "negative"},
                {index =   3, charge =   0.3},
                {index =   5, charge =   0.5},
                {index =   7, charge =   0.7},
                {index = 100, charge = "toomuch"},
            ]
        )"_toml;

        const auto pot = mjolnir::read_debye_huckel_potential<traits_type>(v);

        const auto ignore_within = pot.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(pot.participants().size() ==   6u);
        BOOST_TEST(pot.participants().at(0)  ==   0u);
        BOOST_TEST(pot.participants().at(1)  ==   1u);
        BOOST_TEST(pot.participants().at(2)  ==   3u);
        BOOST_TEST(pot.participants().at(3)  ==   5u);
        BOOST_TEST(pot.participants().at(4)  ==   7u);
        BOOST_TEST(pot.participants().at(5)  == 100u);

        BOOST_TEST(pot.charges().at(  0)  == real_type(  1.0), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(  1)  == real_type( -1.0), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(  3)  == real_type(  0.3), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(  5)  == real_type(  0.5), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(  7)  == real_type(  0.7), tolerance<real_type>());
        BOOST_TEST(pot.charges().at(100)  == real_type(100.0), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_debye_huckel_ignore_self, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_debye_huckel.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction                     = "Pair"
            potential                       = "DebyeHuckel"
            spatial_partition.type          = "Nothing"
            ignore.molecule                 = "Self"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = [
                {index =   0, charge = 1.0},
            ]
        )"_toml;

        const auto pot = mjolnir::read_debye_huckel_potential<traits_type>(v);

        const auto ignore_within = pot.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST( pot.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST( pot.exclusion_list().is_ignored_molecule(1, 1));

        // by default, no group is ignored
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein1", "protein1"));
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein1", "protein2"));
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein2", "protein2"));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_debye_huckel_ignore_others, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_debye_huckel.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction                     = "Pair"
            potential                       = "DebyeHuckel"
            spatial_partition.type          = "Nothing"
            ignore.molecule                 = "Others"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = [
                {index =   0, charge = 1.0},
            ]
        )"_toml;

        const auto pot = mjolnir::read_debye_huckel_potential<traits_type>(v);

        const auto ignore_within = pot.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST( pot.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(1, 1));

        // by default, no group is ignored
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein1", "protein1"));
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein1", "protein2"));
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein2", "protein2"));
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(read_debye_huckel_ignore_group, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_debye_huckel.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction                     = "Pair"
            potential                       = "DebyeHuckel"
            spatial_partition.type          = "Nothing"
            ignore.molecule                 = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            ignore.group.inter      = [
                ["protein1", "protein2"], # between these
                ["protein1", "protein3"],
            ]
            parameters = [
                {index =   0, charge = 1.0},
            ]
        )"_toml;

        const auto pot = mjolnir::read_debye_huckel_potential<traits_type>(v);

        const auto ignore_within = pot.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(1, 1));

        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein1", "protein1"));
        BOOST_TEST( pot.exclusion_list().is_ignored_group("protein1", "protein2"));
        BOOST_TEST( pot.exclusion_list().is_ignored_group("protein1", "protein3"));

        BOOST_TEST( pot.exclusion_list().is_ignored_group("protein2", "protein1"));
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein2", "protein2"));
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein2", "protein3"));

        BOOST_TEST( pot.exclusion_list().is_ignored_group("protein3", "protein1"));
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein3", "protein2"));
        BOOST_TEST(!pot.exclusion_list().is_ignored_group("protein3", "protein3"));

    }
}
