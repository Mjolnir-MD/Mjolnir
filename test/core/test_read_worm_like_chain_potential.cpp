#define BOOST_TEST_MODULE "test_read_worm_like_chain_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_local_potential.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_worm_like_chain_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_worm_like_chain.log");

    using namespace toml::literals;
    using real_type = T;
    const toml::value env;
    const auto param = u8R"(
        indices = [1, 2]
        p       = 3.9
        lc      = 19.0
    )"_toml;

    const auto w = mjolnir::read_worm_like_chain_potential<real_type>(param, env);
    BOOST_TEST(w.p()    == real_type(3.9),  tolerance<real_type>());
    BOOST_TEST(w.lc()   == real_type(19.0), tolerance<real_type>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_worm_like_chain_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_worm_like_chain.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto env = u8R"(
            indices = [1, 2]
            p       = 3.9
            lc      = 19.0
        )"_toml;
        const auto param = u8R"(
            indices = "indices"
            p       = "p"
            lc      = "lc"
        )"_toml;

        const auto w = mjolnir::read_worm_like_chain_potential<real_type>(param, env);
        BOOST_TEST(w.p()     == real_type(3.9),  tolerance<real_type>());
        BOOST_TEST(w.lc()    == real_type(19.0),  tolerance<real_type>());
    }
}
