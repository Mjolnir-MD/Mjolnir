#define BOOST_TEST_MODULE "test_read_flexible_local_angle_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_flexible_local_angle_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_angle.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value env;
        const toml::value v = u8R"(
            indices = [1, 2, 3]
            k       = 3.14
            y       = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
            d2y     = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        )"_toml;

        const auto g = mjolnir::read_flexible_local_angle_potential<real_type>(v, env);
        BOOST_TEST(g.k()      == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(g.y()[0]   == real_type( 1.0), tolerance<real_type>());
        BOOST_TEST(g.y()[1]   == real_type( 2.0), tolerance<real_type>());
        BOOST_TEST(g.y()[2]   == real_type( 3.0), tolerance<real_type>());
        BOOST_TEST(g.y()[3]   == real_type( 4.0), tolerance<real_type>());
        BOOST_TEST(g.y()[4]   == real_type( 5.0), tolerance<real_type>());
        BOOST_TEST(g.y()[5]   == real_type( 6.0), tolerance<real_type>());
        BOOST_TEST(g.y()[6]   == real_type( 7.0), tolerance<real_type>());
        BOOST_TEST(g.y()[7]   == real_type( 8.0), tolerance<real_type>());
        BOOST_TEST(g.y()[8]   == real_type( 9.0), tolerance<real_type>());
        BOOST_TEST(g.y()[9]   == real_type(10.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[0] == real_type( 1.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[1] == real_type( 2.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[2] == real_type( 3.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[3] == real_type( 4.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[4] == real_type( 5.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[5] == real_type( 6.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[6] == real_type( 7.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[7] == real_type( 8.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[8] == real_type( 9.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[9] == real_type(10.0), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_flexible_local_angle_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_angle.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value env = u8R"(
            indices = [1, 2, 3]
            k       = 3.14
            y       = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
            d2y     = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        )"_toml;
        const toml::value v = u8R"(
            indices = "indices"
            k       = "k"
            y       = "y"
            d2y     = "d2y"
        )"_toml;

        const auto g = mjolnir::read_flexible_local_angle_potential<real_type>(v, env);
        BOOST_TEST(g.k()      == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(g.y()[0]   == real_type( 1.0), tolerance<real_type>());
        BOOST_TEST(g.y()[1]   == real_type( 2.0), tolerance<real_type>());
        BOOST_TEST(g.y()[2]   == real_type( 3.0), tolerance<real_type>());
        BOOST_TEST(g.y()[3]   == real_type( 4.0), tolerance<real_type>());
        BOOST_TEST(g.y()[4]   == real_type( 5.0), tolerance<real_type>());
        BOOST_TEST(g.y()[5]   == real_type( 6.0), tolerance<real_type>());
        BOOST_TEST(g.y()[6]   == real_type( 7.0), tolerance<real_type>());
        BOOST_TEST(g.y()[7]   == real_type( 8.0), tolerance<real_type>());
        BOOST_TEST(g.y()[8]   == real_type( 9.0), tolerance<real_type>());
        BOOST_TEST(g.y()[9]   == real_type(10.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[0] == real_type( 1.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[1] == real_type( 2.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[2] == real_type( 3.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[3] == real_type( 4.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[4] == real_type( 5.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[5] == real_type( 6.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[6] == real_type( 7.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[7] == real_type( 8.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[8] == real_type( 9.0), tolerance<real_type>());
        BOOST_TEST(g.d2y()[9] == real_type(10.0), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_local_potential_flexible_local_angle_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_angle.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            parameters = [
                {indices = [1, 2, 3], k = 3.14, y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], d2y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<3,
              mjolnir::FlexibleLocalAnglePotential<real_type>>(v);

        const std::array<std::size_t, 3> ref_idx{{1, 2, 3}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()      == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[0]   == real_type( 1.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[1]   == real_type( 2.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[2]   == real_type( 3.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[3]   == real_type( 4.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[4]   == real_type( 5.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[5]   == real_type( 6.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[6]   == real_type( 7.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[7]   == real_type( 8.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[8]   == real_type( 9.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[9]   == real_type(10.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[0] == real_type( 1.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[1] == real_type( 2.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[2] == real_type( 3.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[3] == real_type( 4.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[4] == real_type( 5.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[5] == real_type( 6.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[6] == real_type( 7.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[7] == real_type( 8.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[8] == real_type( 9.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[9] == real_type(10.0), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_local_potential_flexible_local_angle_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_angle.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            env.pi  = 3.14
            env.y   = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
            env.d2y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
            parameters = [
                {indices = [1, 2, 3], k = "pi", y = "y", d2y = "d2y"}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<3,
              mjolnir::FlexibleLocalAnglePotential<real_type>>(v);

        const std::array<std::size_t, 3> ref_idx{{1, 2, 3}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()      == real_type(3.14), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[0]   == real_type( 1.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[1]   == real_type( 2.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[2]   == real_type( 3.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[3]   == real_type( 4.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[4]   == real_type( 5.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[5]   == real_type( 6.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[6]   == real_type( 7.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[7]   == real_type( 8.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[8]   == real_type( 9.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.y()[9]   == real_type(10.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[0] == real_type( 1.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[1] == real_type( 2.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[2] == real_type( 3.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[3] == real_type( 4.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[4] == real_type( 5.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[5] == real_type( 6.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[6] == real_type( 7.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[7] == real_type( 8.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[8] == real_type( 9.0), tolerance<real_type>());
        BOOST_TEST(g.at(0).second.d2y()[9] == real_type(10.0), tolerance<real_type>());
    }
}
