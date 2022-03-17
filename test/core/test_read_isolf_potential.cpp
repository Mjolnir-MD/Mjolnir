#define BOOST_TEST_MODULE "test_read_isolf_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_isolf_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_isolf.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "iSoLFAttractive"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond  = 1
            ignore.particles_within.angle = 1
            env.large = 100.0
            parameters = [
                {index =   0, sigma =     2.0, epsilon = 1.5, omega = 1.0},
                {index =   1, sigma =     2.0, epsilon = 1.5, omega = 1.0},
                {index =   2, sigma =     5.0, epsilon = 0.5, omega = 1.0},
                {index =   7, sigma =     7.0, epsilon = 0.7, omega = 1.0},
                {index = 100, sigma = "large", epsilon = 0.1, omega = 1.0},
            ]
        )"_toml;

        const auto pot_para = mjolnir::read_isolf_potential<traits_type>(v);
        const auto& para = dynamic_cast<
            mjolnir::iSoLFAttractiveParameterList<traits_type> const&
            >(pot_para.second.cref());

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")  == 1ul);
        BOOST_TEST(within.at("angle") == 1ul);

        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(1, 1));

        BOOST_TEST(para.participants().size() ==   5u);
        BOOST_TEST(para.participants().at(0)  ==   0u);
        BOOST_TEST(para.participants().at(1)  ==   1u);
        BOOST_TEST(para.participants().at(2)  ==   2u);
        BOOST_TEST(para.participants().at(3)  ==   7u);
        BOOST_TEST(para.participants().at(4)  == 100u);

        BOOST_TEST(para.parameters().at(  0).sigma == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  1).sigma == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  2).sigma == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  7).sigma == real_type(  7.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(100).sigma == real_type(100.0), tolerance<real_type>());

        BOOST_TEST(para.parameters().at(  0).epsilon == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  1).epsilon == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  2).epsilon == real_type(0.5), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  7).epsilon == real_type(0.7), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(100).epsilon == real_type(0.1), tolerance<real_type>());

        BOOST_TEST(para.parameters().at(  0).omega == real_type(1.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  1).omega == real_type(1.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  2).omega == real_type(1.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  7).omega == real_type(1.0), tolerance<real_type>());
        BOOST_TEST(para.parameters().at(100).omega == real_type(1.0), tolerance<real_type>());
    }
}
