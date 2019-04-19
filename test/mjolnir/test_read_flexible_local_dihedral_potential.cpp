#define BOOST_TEST_MODULE "test_read_flexible_local_dihedral_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_local_potential.hpp>

BOOST_AUTO_TEST_CASE(read_flexible_local_dihedral_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_dihedral.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            indices = [1, 2, 3, 4]
            k       = 3.14
            coef    = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        )"_toml;

        const auto g = mjolnir::read_flexible_local_dihedral_potential<real_type>(v);
        BOOST_TEST(g.k()       == 3.14, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[0] ==  1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[1] ==  2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[2] ==  3.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[3] ==  4.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[4] ==  5.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[5] ==  6.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[6] ==  7.0, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_flexible_local_dihedral_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_dihedral.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            indices = [1, 2, 3, 4]
            k       = 3.14
            coef    = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        )"_toml;

        const auto g = mjolnir::read_flexible_local_dihedral_potential<real_type>(v);
        BOOST_TEST(g.k()     == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[0] ==  1.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[1] ==  2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[2] ==  3.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[3] ==  4.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[4] ==  5.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[5] ==  6.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.coef()[6] ==  7.0f, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_local_potential_flexible_local_dihedral_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_dihedral.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            parameters = [
                {indices = [1, 2, 3, 4], k = 3.14, coef = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<4,
              mjolnir::FlexibleLocalDihedralPotential<real_type>>(v);

        const std::array<std::size_t, 4> ref_idx{{1, 2, 3, 4}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()       == 3.14, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[0] ==  1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[1] ==  2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[2] ==  3.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[3] ==  4.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[4] ==  5.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[5] ==  6.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[6] ==  7.0, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_local_potential_flexible_local_dihedral_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_flexible_local_dihedral.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            parameters = [
                {indices = [1, 2, 3, 4], k = 3.14, coef = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<4,
              mjolnir::FlexibleLocalDihedralPotential<real_type>>(v);

        const std::array<std::size_t, 4> ref_idx{{1, 2, 3, 4}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k()       == 3.14f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[0] ==  1.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[1] ==  2.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[2] ==  3.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[3] ==  4.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[4] ==  5.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[5] ==  6.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.coef()[6] ==  7.0f, boost::test_tools::tolerance(tol));
    }
}
