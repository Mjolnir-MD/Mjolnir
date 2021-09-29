#define BOOST_TEST_MODULE "test_read_hard_core_excluded_volume_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_hard_core_excluded_volume_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_hard_core_excluded_volume.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            epsilon                 = 3.14
            cutoff                  = 2.71
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters  = [
                {index = 0, core_radius = 2.0, soft_shell_thickness = 3.0},
                {index = 1, core_radius = 2.0, soft_shell_thickness = 2.0},
                {index = 3, core_radius = 3.0, soft_shell_thickness = 100.0},
            ]
        )"_toml;

        const auto pot_para = mjolnir::read_hard_core_excluded_volume_potential<traits_type>(v);
        const auto& pot  = pot_para.first;
        const auto& para = dynamic_cast<mjolnir::HardCoreExcludedVolumeParameterList<traits_type> const&>(pot_para.second.cref());

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(para.participants().size() == 3u);
        BOOST_TEST(para.participants().at(0)  == 0u);
        BOOST_TEST(para.participants().at(1)  == 1u);
        BOOST_TEST(para.participants().at(2)  == 3u);

        // soft_shell_thickness (sigma)
        BOOST_TEST(para.parameters().at(0).sigma == real_type(  3.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(1).sigma == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(3).sigma == real_type(100.0), tolerance<real_type>());

        // core_radius (r0)
        BOOST_TEST(para.parameters().at(0).radius == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(1).radius == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(3).radius == real_type(3.0), tolerance<real_type>());

        BOOST_TEST(pot.epsilon()      == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(pot.cutoff_ratio() == real_type(2.71), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_hard_core_excluded_volume_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_hard_core_excluded_volume.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            epsilon                 = 3.14
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            env.two      = 2.0
            env.three    = 3.0
            env.toolarge = 100.0
            parameters  = [
                {index =   0, core_radius =   "two", soft_shell_thickness = 3.0},
                {index =   1, core_radius =     2.0, soft_shell_thickness = "two"},
                {index =   3, core_radius = "three", soft_shell_thickness = "toolarge"},
            ]
        )"_toml;

        const auto pot_para = mjolnir::read_hard_core_excluded_volume_potential<traits_type>(v);
        const auto& pot  = pot_para.first;
        const auto& para = dynamic_cast<mjolnir::HardCoreExcludedVolumeParameterList<traits_type> const&>(pot_para.second.cref());

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(para.participants().size() == 3u);
        BOOST_TEST(para.participants().at(0)  == 0u);
        BOOST_TEST(para.participants().at(1)  == 1u);
        BOOST_TEST(para.participants().at(2)  == 3u);

        // soft_shell_thickness (sigma)
        BOOST_TEST(para.parameters().at(0).sigma == real_type(  3.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(1).sigma == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(3).sigma == real_type(100.0), tolerance<real_type>());

        // core radius (r0)
        BOOST_TEST(para.parameters().at(0).radius == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(1).radius == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(3).radius == real_type(3.0), tolerance<real_type>());

        BOOST_TEST(pot.epsilon()      == real_type(3.14), tolerance<real_type>());
    }
}
