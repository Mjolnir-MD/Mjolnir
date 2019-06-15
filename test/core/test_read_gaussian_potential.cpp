#define BOOST_TEST_MODULE "test_read_gaussian_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_gaussian_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_gaussian.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value env;
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
            sigma   = 0.577
            v0      = 2.71
        )"_toml;

        const auto g = mjolnir::read_gaussian_potential<real_type>(v, env);
        BOOST_TEST(g.k()     == real_type(3.14),  tolerance<real_type>());
        BOOST_TEST(g.sigma() == real_type(0.577), tolerance<real_type>());
        BOOST_TEST(g.v0()    == real_type(2.71),  tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_gaussian_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_gaussian.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto env = u8R"(
            indices = [1, 2]
            k       = 3.14
            sigma   = 0.577
            v0      = 2.71
        )"_toml;
        const auto v = u8R"(
            indices = "indices"
            k       = "k"
            sigma   = "sigma"
            v0      = "v0"
        )"_toml;

        const auto g = mjolnir::read_gaussian_potential<real_type>(v, env);
        BOOST_TEST(g.k()     == real_type(3.14),  tolerance<real_type>());
        BOOST_TEST(g.sigma() == real_type(0.577), tolerance<real_type>());
        BOOST_TEST(g.v0()    == real_type(2.71),  tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_local_potential_gaussian_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_gaussian.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            parameters = [
                {indices = [1,2], k = 3.14, sigma = 0.577, v0 = 2.71}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
            mjolnir::GaussianPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()     == real_type(3.14),  tolerance<real_type>());
        BOOST_TEST(g.at(0).second.sigma() == real_type(0.577), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.v0()    == real_type(2.71),  tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_local_potential_gaussian_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_gaussian.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            env.indices = [1, 2]
            env.pi      = 3.14
            env.e       = 2.71
            env.gamma   = 0.577
            parameters = [
                {indices = [1,2], k = "pi", sigma = "gamma", v0 = "e"}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
            mjolnir::GaussianPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()     == real_type(3.14),  tolerance<real_type>());
        BOOST_TEST(g.at(0).second.sigma() == real_type(0.577), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.v0()    == real_type(2.71),  tolerance<real_type>());
    }
}
