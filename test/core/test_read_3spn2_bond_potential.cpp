#define BOOST_TEST_MODULE "test_read_3spn2_bond_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_3spn2_bond_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_3spn2_bond.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value env;
        const auto v = u8R"(
            indices = [1, 2]
            k       = 3.14
            v0      = 2.71
        )"_toml;

        const auto pot = mjolnir::read_3spn2_bond_potential<real_type>(v, env);
        BOOST_TEST(pot.k()  == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(pot.v0() == real_type(2.71), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_3spn2_bond_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_3spn2_bond.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto env = u8R"(
            indices = [1, 2]
            k       = 3.14
            v0      = 2.71
        )"_toml;
        const auto v = u8R"(
            indices = "indices"
            k       = "k"
            v0      = "v0"
        )"_toml;
        const auto pot = mjolnir::read_3spn2_bond_potential<real_type>(v, env);
        BOOST_TEST(pot.k()  == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(pot.v0() == real_type(2.71), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_local_potential_3spn2_bond_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_3spn2_bond.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            parameters = [
                {indices = [1,2], k = 3.14, v0 = 2.71}
            ]
        )"_toml;

        const auto pot = mjolnir::read_local_potential<2,
            mjolnir::ThreeSPN2BondPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(pot.size() == 1u);
        BOOST_TEST(pot.at(0).first == ref_idx);
        BOOST_TEST(pot.at(0).second.k()  == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(pot.at(0).second.v0() == real_type(2.71), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_local_potential_3spn2_bond_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_3spn2_bond.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            env.indices = [1, 2]
            env.pi      = 3.14
            env.e       = 2.71
            parameters = [
                {indices = "indices", k = "pi", v0 = "e"}
            ]
        )"_toml;

        const auto pot = mjolnir::read_local_potential<2,
            mjolnir::ThreeSPN2BondPotential<real_type>>(v);

        const std::array<std::size_t, 2> ref_idx{{1, 2}};

        BOOST_TEST(pot.size() == 1u);
        BOOST_TEST(pot.at(0).first == ref_idx);
        BOOST_TEST(pot.at(0).second.k()  == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(pot.at(0).second.v0() == real_type(2.71), tolerance<real_type>());
    }
}
