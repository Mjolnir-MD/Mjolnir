#define BOOST_TEST_MODULE "test_read_tabulated_wca_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_global_potential.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_tabulated_wca_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_tabulated_wca.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "WCA"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            table.A.A = {sigma =   2.0, epsilon = 1.5}
            table.A.B = {sigma =   5.0, epsilon = 0.5}
            table.B.B = {sigma = 100.0, epsilon = 0.1}
            parameters = [
                {index =   0, name = "A"},
                {index =   1, name = "B"},
                {index =   2, name = "A"},
                {index =   3, name = "B"},
                {index =   5, name = "A"},
                {index =   7, name = "B"},
                {index = 100, name = "A"},
            ]
        )"_toml;

        const auto pot = mjolnir::read_tabulated_wca_potential<traits_type>(v);

        const auto ignore_within = pot.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(1, 1));

        BOOST_TEST(pot.participants().size() ==   7u);
        BOOST_TEST(pot.participants().at(0)  ==   0u);
        BOOST_TEST(pot.participants().at(1)  ==   1u);
        BOOST_TEST(pot.participants().at(2)  ==   2u);
        BOOST_TEST(pot.participants().at(3)  ==   3u);
        BOOST_TEST(pot.participants().at(4)  ==   5u);
        BOOST_TEST(pot.participants().at(5)  ==   7u);
        BOOST_TEST(pot.participants().at(6)  == 100u);

        BOOST_TEST(pot.parameters().at(  0)  == "A");
        BOOST_TEST(pot.parameters().at(  1)  == "B");
        BOOST_TEST(pot.parameters().at(  2)  == "A");
        BOOST_TEST(pot.parameters().at(  3)  == "B");
        BOOST_TEST(pot.parameters().at(  5)  == "A");
        BOOST_TEST(pot.parameters().at(  7)  == "B");
        BOOST_TEST(pot.parameters().at(100)  == "A");

        BOOST_TEST(pot.prepare_params(0, 1).first == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(0, 2).first == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(1, 2).first == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(1, 3).first == real_type(100.0), tolerance<real_type>());

        BOOST_TEST(pot.prepare_params(0, 1).second == real_type(  0.5), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(0, 2).second == real_type(  1.5), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(1, 2).second == real_type(  0.5), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(1, 3).second == real_type(  0.1), tolerance<real_type>());

    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_tabulated_wca_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_tabulated_wca.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "WCA"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1

            env.100   = 100.0
            env.small = 0.1
            env.half  = 0.5
            table.A.A = {sigma =   2.0, epsilon = 1.5}
            table.A.B = {sigma =   5.0, epsilon = "half"}
            table.B.B = {sigma = "100", epsilon = "small"}
            parameters = [
                {index =   0, name = "A"},
                {index =   1, name = "B"},
                {index =   2, name = "A"},
                {index =   3, name = "B"},
                {index =   5, name = "A"},
                {index =   7, name = "B"},
                {index = 100, name = "A"},
            ]
        )"_toml;

        const auto pot = mjolnir::read_tabulated_wca_potential<traits_type>(v);

        const auto ignore_within = pot.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(1, 1));

        BOOST_TEST(pot.participants().size() ==   7u);
        BOOST_TEST(pot.participants().at(0)  ==   0u);
        BOOST_TEST(pot.participants().at(1)  ==   1u);
        BOOST_TEST(pot.participants().at(2)  ==   2u);
        BOOST_TEST(pot.participants().at(3)  ==   3u);
        BOOST_TEST(pot.participants().at(4)  ==   5u);
        BOOST_TEST(pot.participants().at(5)  ==   7u);
        BOOST_TEST(pot.participants().at(6)  == 100u);

        BOOST_TEST(pot.parameters().at(  0)  == "A");
        BOOST_TEST(pot.parameters().at(  1)  == "B");
        BOOST_TEST(pot.parameters().at(  2)  == "A");
        BOOST_TEST(pot.parameters().at(  3)  == "B");
        BOOST_TEST(pot.parameters().at(  5)  == "A");
        BOOST_TEST(pot.parameters().at(  7)  == "B");
        BOOST_TEST(pot.parameters().at(100)  == "A");

        BOOST_TEST(pot.prepare_params(0, 1).first == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(0, 2).first == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(1, 2).first == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(1, 3).first == real_type(100.0), tolerance<real_type>());

        BOOST_TEST(pot.prepare_params(0, 1).second == real_type(  0.5), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(0, 2).second == real_type(  1.5), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(1, 2).second == real_type(  0.5), tolerance<real_type>());
        BOOST_TEST(pot.prepare_params(1, 3).second == real_type(  0.1), tolerance<real_type>());


    }
}
