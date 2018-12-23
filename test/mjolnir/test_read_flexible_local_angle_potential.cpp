#define BOOST_TEST_MODULE "test_read_flexible_local_angle_potential"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/input/read_potential.hpp>

BOOST_AUTO_TEST_CASE(read_flexible_local_angle_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_angle.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",   toml::value(3.14)},
            {"y",   toml::value{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}},
            {"d2y", toml::value{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}},
        };
        const auto g = mjolnir::read_flexible_local_angle_potential<real_type>(v);
        BOOST_TEST(g.k()   == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[0] ==  1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[1] ==  2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[2] ==  3.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[3] ==  4.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[4] ==  5.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[5] ==  6.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[6] ==  7.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[7] ==  8.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[8] ==  9.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[9] == 10.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[0] ==  1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[1] ==  2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[2] ==  3.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[3] ==  4.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[4] ==  5.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[5] ==  6.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[6] ==  7.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[7] ==  8.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[8] ==  9.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[9] == 10.0, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_flexible_local_angle_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_angle.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        const toml::value v = toml::table{
            {"indices", toml::value({1, 2})},
            {"k",   toml::value(3.14)},
            {"y",   toml::value{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}},
            {"d2y", toml::value{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}},
        };
        const auto g = mjolnir::read_flexible_local_angle_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[0] ==  1.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[1] ==  2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[2] ==  3.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[3] ==  4.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[4] ==  5.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[5] ==  6.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[6] ==  7.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[7] ==  8.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[8] ==  9.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.y()[9] == 10.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[0] ==  1.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[1] ==  2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[2] ==  3.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[3] ==  4.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[4] ==  5.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[5] ==  6.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[6] ==  7.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[7] ==  8.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[8] ==  9.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.d2y()[9] == 10.0f, boost::test_tools::tolerance(tol));    }
}
