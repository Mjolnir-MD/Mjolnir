#define BOOST_TEST_MODULE "test_read_clementi_dihedral_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_local_potential.hpp>

BOOST_AUTO_TEST_CASE(read_clementi_dihedral_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_clementi_dihedral.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            indices = [1, 2, 3, 4]
            k1 = 3.14
            k3 = 0.577
            v0 = 2.71
        )"_toml;

        const auto g = mjolnir::read_clementi_dihedral_potential<real_type>(v);
        BOOST_TEST(g.k1() == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.k3() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0() == 2.71,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_clementi_dihedral_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_clementi_dihedral.log");
    using real_type = float;
    constexpr real_type tol = 1e-4;

    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            indices = [1, 2, 3, 4]
            k1 = 3.14
            k3 = 0.577
            v0 = 2.71
        )"_toml;

        const auto g = mjolnir::read_clementi_dihedral_potential<real_type>(v);
        BOOST_TEST(g.k1() == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.k3() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.v0() == 2.71f,  boost::test_tools::tolerance(tol));
    }
}

// ---------------------------------------------------------------------------
// read_local_potential

BOOST_AUTO_TEST_CASE(read_local_potential_clementi_dihedral_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_clementi_dihedral.log");

    using real_type = double;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            parameters = [
                {indices = [1, 2, 3, 4], k1 = 3.14, k3 = 0.577, v0 = 2.71}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<4,
              mjolnir::ClementiDihedralPotential<real_type>>(v);

        const std::array<std::size_t, 4> ref_idx{{1, 2, 3, 4}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k1() == 3.14,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.k3() == 0.577, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0() == 2.71,  boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_local_potential_clementi_dihedral_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_clementi_dihedral.log");

    using real_type = float;
    constexpr real_type tol = 1e-4;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            parameters = [
                {indices = [1, 2, 3, 4], k1 = 3.14, k3 = 0.577, v0 = 2.71}
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<4,
              mjolnir::ClementiDihedralPotential<real_type>>(v);

        const std::array<std::size_t, 4> ref_idx{{1, 2, 3, 4}};

        BOOST_TEST(g.size() == 1u);
        BOOST_TEST(g.at(0).first == ref_idx);
        BOOST_TEST(g.at(0).second.k1() == 3.14f,  boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.k3() == 0.577f, boost::test_tools::tolerance(tol));
        BOOST_TEST(g.at(0).second.v0() == 2.71f,  boost::test_tools::tolerance(tol));
    }
}
