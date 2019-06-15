#define BOOST_TEST_MODULE "test_read_harmonic_restraint_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_external_potential.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_harmonic_restraint_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic_restraint.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "PositionRestraint"
            potential   = "Harmonic"
            parameters  = [
                {position = [0.0, 0.0, 0.0], index = 0, k = 10.0, v0 = 0.0},
                {position = [1.0, 2.0, 2.0], index = 1, k = 20.0, v0 = 1.0},
            ]
        )"_toml;

        const auto h = mjolnir::read_harmonic_restraint_potential<real_type>(v);

        BOOST_TEST(h.parameters().size() == 2u);
        BOOST_TEST(h.parameters().at(0).first  == real_type(10.0), tolerance<real_type>());
        BOOST_TEST(h.parameters().at(0).second == real_type( 0.0), tolerance<real_type>());
        BOOST_TEST(h.parameters().at(1).first  == real_type(20.0), tolerance<real_type>());
        BOOST_TEST(h.parameters().at(1).second == real_type( 1.0), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_harmonic_restraint_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic_restraint.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "PositionRestraint"
            potential   = "Harmonic"
            env.strong = 100.0
            env.weak   =   0.1
            env.just   =   0.0
            env.far    = 100.0
            parameters  = [
                {position = [0.0, 0.0, 0.0], index = 0, k = "strong", v0 = "just"},
                {position = [1.0, 2.0, 2.0], index = 1, k = "weak",   v0 = "far"},
            ]
        )"_toml;

        const auto h = mjolnir::read_harmonic_restraint_potential<real_type>(v);

        BOOST_TEST(h.parameters().size() == 2u);
        BOOST_TEST(h.parameters().at(0).first  == real_type(100.0), tolerance<real_type>());
        BOOST_TEST(h.parameters().at(0).second == real_type(  0.0), tolerance<real_type>());
        BOOST_TEST(h.parameters().at(1).first  == real_type(  0.1), tolerance<real_type>());
        BOOST_TEST(h.parameters().at(1).second == real_type(100.0), tolerance<real_type>());
    }
}
