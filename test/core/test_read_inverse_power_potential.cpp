#define BOOST_TEST_MODULE "test_read_inverse_power_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_inverse_power_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_inverse_power.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::InversePowerPotential<real_type>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            epsilon                 = 3.14
            n                       = 5
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters  = [
                {index =   0, radius =   2.0},
                {index =   1, radius =   2.0},
                {index =   3, radius =   3.0},
                {index =   5, radius =   5.0},
                {index =   7, radius =   7.0},
                {index = 100, radius = 100.0},
            ]
        )"_toml;

        const auto pot_para = mjolnir::read_inverse_power_potential<traits_type>(v);
        const auto& pot  = pot_para.first;
        const auto& para = dynamic_cast<mjolnir::InversePowerParameterList<traits_type> const&>(pot_para.second.cref());

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")     == 3ul);
        BOOST_TEST(within.at("contact")  == 1ul);

        BOOST_TEST(para.participants().size() ==   6u);
        BOOST_TEST(para.participants().at(0)  ==   0u);
        BOOST_TEST(para.participants().at(1)  ==   1u);
        BOOST_TEST(para.participants().at(2)  ==   3u);
        BOOST_TEST(para.participants().at(3)  ==   5u);
        BOOST_TEST(para.participants().at(4)  ==   7u);
        BOOST_TEST(para.participants().at(5)  == 100u);

        BOOST_TEST(para.parameters().at(  0).sigma == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  1).sigma == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  3).sigma == real_type(  3.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  5).sigma == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  7).sigma == real_type(  7.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(100).sigma == real_type(100.0), tolerance<real_type>());


        BOOST_TEST(pot.epsilon() == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(pot.n()       == 5);
        BOOST_TEST(pot.cutoff_ratio() == potential_type::default_cutoff(5), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_inverse_power_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_inverse_power.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::InversePowerPotential<real_type>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            epsilon                 = 3.14
            n                       = 5
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            env.five     = 5.0
            env.seven    = 7.0
            env.toolarge = 100.0
            parameters  = [
                {index =   0, radius =   2.0},
                {index =   1, radius =   2.0},
                {index =   3, radius =   3.0},
                {index =   5, radius =   "five"},
                {index =   7, radius =   "seven"},
                {index = 100, radius =   "toolarge"},
            ]
        )"_toml;

        const auto pot_para = mjolnir::read_inverse_power_potential<traits_type>(v);
        const auto& pot  = pot_para.first;
        const auto& para = dynamic_cast<mjolnir::InversePowerParameterList<traits_type> const&>(pot_para.second.cref());

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")     == 3ul);
        BOOST_TEST(within.at("contact")  == 1ul);

        BOOST_TEST(para.participants().size() ==   6u);
        BOOST_TEST(para.participants().at(0)  ==   0u);
        BOOST_TEST(para.participants().at(1)  ==   1u);
        BOOST_TEST(para.participants().at(2)  ==   3u);
        BOOST_TEST(para.participants().at(3)  ==   5u);
        BOOST_TEST(para.participants().at(4)  ==   7u);
        BOOST_TEST(para.participants().at(5)  == 100u);

        BOOST_TEST(para.parameters().at(  0).sigma == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  1).sigma == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  3).sigma == real_type(  3.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  5).sigma == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  7).sigma == real_type(  7.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(100).sigma == real_type(100.0), tolerance<real_type>());


        BOOST_TEST(pot.epsilon() == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(pot.n()       == 5);
        BOOST_TEST(pot.cutoff_ratio() == potential_type::default_cutoff(5), tolerance<real_type>());
    }
}
