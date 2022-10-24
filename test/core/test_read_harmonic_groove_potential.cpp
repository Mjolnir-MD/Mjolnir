#define BOOST_TEST_MODULE "test_read_harmonic_groove_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_external_potential.hpp>
#include <tuple>
#include <limits>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e08;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_harmonic_groove_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic_groove.log");

    using real_type = T;
    {
        using namespace toml::literals;
        // TODO: enable to remove shape.margin line
        const toml::value v = u8R"(
            interaction    = "Distance"
            potential      = "HarmonicGroove"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "+X"
            shape.position = 1.0
            shape.margin   = 1.0
            parameters     = [
                {index =  0, k = 2.0, v0 = 1.5},
                {index =  9, k = 2.0, v0 = 1.5},
            ]
        )"_toml;

        const auto g = mjolnir::read_harmonic_groove_potential<real_type>(v);

        BOOST_TEST(g.participants().size() ==  2u);
        BOOST_TEST(g.participants().at(0)  ==  0u);
        BOOST_TEST(g.participants().at(1)  ==  9u);

        BOOST_TEST(g.parameters().at( 0).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 9).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 0).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 9).second == real_type(1.5), tolerance<real_type>());
    }

    {
        using namespace toml::literals;
        // TODO: enable to remove shape.margin line
        const toml::value v = u8R"(
            interaction    = "Distance"
            potential      = "HarmonicGroove"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "+X"
            shape.position = 1.0
            shape.margin   = 1.0
            parameters     = [
                {index =  0, k = 2.0, v0 = 1.5},
                {index =  9, k = 2.0, v0 = 1.5},
            ]
        )"_toml;

        const auto g = mjolnir::read_harmonic_groove_potential<real_type>(v);

        BOOST_TEST(g.participants().size() ==  2u);
        BOOST_TEST(g.participants().at(0)  ==  0u);
        BOOST_TEST(g.participants().at(1)  ==  9u);

        BOOST_TEST(g.parameters().at( 0).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 9).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 0).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 9).second == real_type(1.5), tolerance<real_type>());
    }

}
