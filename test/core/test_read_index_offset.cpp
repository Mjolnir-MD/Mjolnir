#define BOOST_TEST_MODULE "test_read_harmonic_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_local_potential.hpp>
#include <mjolnir/input/read_global_potential.hpp>
#include <mjolnir/input/read_external_potential.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_index_offset_local_potential, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_index_offset.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            parameters = [
                {indices = [1, 2],              k = 3.14, v0 = 2.71},
                {indices = [1, 2], offset =  1, k = 3.14, v0 = 2.71},
                {indices = [1, 2], offset = -1, k = 3.14, v0 = 2.71},
            ]
        )"_toml;

        const auto g = mjolnir::read_local_potential<2,
            mjolnir::HarmonicPotential<real_type>>(v);

        BOOST_TEST(g.size() == 3u);
        BOOST_TEST(g.at(0).first[0] == 1u);
        BOOST_TEST(g.at(0).first[1] == 2u);

        BOOST_TEST(g.at(1).first[0] == 2u);
        BOOST_TEST(g.at(1).first[1] == 3u);

        BOOST_TEST(g.at(2).first[0] == 0u);
        BOOST_TEST(g.at(2).first[1] == 1u);

        BOOST_TEST(g.at(0).second.k()  == real_type(3.14));
        BOOST_TEST(g.at(0).second.v0() == real_type(2.71));
        BOOST_TEST(g.at(1).second.k()  == real_type(3.14));
        BOOST_TEST(g.at(1).second.v0() == real_type(2.71));
        BOOST_TEST(g.at(2).second.k()  == real_type(3.14));
        BOOST_TEST(g.at(2).second.v0() == real_type(2.71));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_index_offset_global_potential, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_index_offset.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "LennardJones"
            spatial_partition.type  = "Naive"
            parameters = [
                {index =   0,               sigma =   2.0, epsilon = 1.5},
                {index =   5, offset =   1, sigma =   5.0, epsilon = 0.5},
                {index =  10, offset =  -1, sigma =   7.0, epsilon = 0.7},
                {index = 100, offset = 100, sigma = 100.0, epsilon = 0.1},
            ]
        )"_toml;

        const auto pot = mjolnir::read_lennard_jones_potential<traits_type>(v);

        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!pot.exclusion_list().is_ignored_molecule(1, 1));

        BOOST_TEST(pot.participants().size() ==   4u);
        BOOST_TEST(pot.participants().at(0)  ==   0u);
        BOOST_TEST(pot.participants().at(1)  ==   6u);
        BOOST_TEST(pot.participants().at(2)  ==   9u);
        BOOST_TEST(pot.participants().at(3)  == 200u);

        BOOST_TEST(pot.parameters().at(  0).first  == real_type(  2.0), tolerance<real_type>());
        BOOST_TEST(pot.parameters().at(  6).first  == real_type(  5.0), tolerance<real_type>());
        BOOST_TEST(pot.parameters().at(  9).first  == real_type(  7.0), tolerance<real_type>());
        BOOST_TEST(pot.parameters().at(200).first  == real_type(100.0), tolerance<real_type>());

        BOOST_TEST(pot.parameters().at(  0).second == real_type(1.5), tolerance<real_type>());
        BOOST_TEST(pot.parameters().at(  6).second == real_type(0.5), tolerance<real_type>());
        BOOST_TEST(pot.parameters().at(  9).second == real_type(0.7), tolerance<real_type>());
        BOOST_TEST(pot.parameters().at(200).second == real_type(0.1), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_index_offset_external_potential, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_index_offset.log");

    using real_type = T;
    using real_type = T;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction    = "Distance"
            potential      = "LennardJonesWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "+X"
            shape.position = 1.0
            shape.margin   = 0.5
            cutoff         = 3.14
            parameters     = [
                {index =  0,               sigma = 1.0, epsilon = 1.0},
                {index =  2, offset =   1, sigma = 2.0, epsilon = 2.0},
                {index =  9, offset =  -1, sigma = 3.0, epsilon = 3.0},
                {index = 10, offset =  10, sigma = 4.0, epsilon = 4.0},
            ]
        )"_toml;

        const auto g = mjolnir::read_lennard_jones_wall_potential<real_type>(v);

        BOOST_TEST(g.participants().size() ==  4u);
        BOOST_TEST(g.participants().at(0)  ==  0u);
        BOOST_TEST(g.participants().at(1)  ==  3u);
        BOOST_TEST(g.participants().at(2)  ==  8u);
        BOOST_TEST(g.participants().at(3)  == 20u);

        BOOST_TEST(g.parameters().at( 0).first  == real_type(1.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 3).first  == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 8).first  == real_type(3.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(20).first  == real_type(4.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 0).second == real_type(1.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 3).second == real_type(2.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at( 8).second == real_type(3.0), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(20).second == real_type(4.0), tolerance<real_type>());
        BOOST_TEST(g.cutoff() == real_type(3.14), tolerance<real_type>());
    }
}
