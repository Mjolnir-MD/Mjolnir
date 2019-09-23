#define BOOST_TEST_MODULE "test_read_implicit_membrane_potential"

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

BOOST_AUTO_TEST_CASE_TEMPLATE(read_implicit_membrane_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_implicit_membrane.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction    = "Distance"
            potential      = "ExcludedVolumeWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "-X"
            shape.position = 1.0
            shape.margin   = 0.5
            thickness      = 3.14
            interaction_magnitude = 6.28
            bend                  = 9.42
            parameters = [
                {index = 0, hydrophobicity = 2.0},
                {index = 1, hydrophobicity = 3.0},
                {index = 5, hydrophobicity = 4.0},
            ]
        )"_toml;

        const auto g = mjolnir::read_implicit_membrane_potential<real_type>(v);

        BOOST_TEST(g.participants().size() == 3u);
        BOOST_TEST(g.participants().at(0)  == 0u);
        BOOST_TEST(g.participants().at(1)  == 1u);
        BOOST_TEST(g.participants().at(2)  == 5u);
        BOOST_TEST(g.parameters().at(0)    == real_type(2.0       ), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(1)    == real_type(3.0       ), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(5)    == real_type(4.0       ), tolerance<real_type>());
        BOOST_TEST(g.half_thick()          == real_type(3.14 * 0.5), tolerance<real_type>());
        BOOST_TEST(g.k()                   == real_type(6.28      ), tolerance<real_type>());
        BOOST_TEST(g.bend()                == real_type(9.42      ), tolerance<real_type>());
        BOOST_TEST(g.cutoff()              == real_type(4.0       ), tolerance<real_type>());
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_implicit_membrane_env, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_implicit_membrane.log");

    using real_type = T;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction    = "Distance"
            potential      = "ExcludedVolumeWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "-X"
            shape.position = 1.0
            shape.margin   = 0.5
            thickness      = 3.14
            cutoff         = 5.0
            interaction_magnitude = 6.28
            bend                  = 9.42
            env.hydrophilic =   0.0
            env.hydrophobic = 100.0
            parameters = [
                {index = 0, hydrophobicity = "hydrophilic"},
                {index = 1, hydrophobicity = "hydrophobic"},
                {index = 5, hydrophobicity = "hydrophobic"},
            ]
        )"_toml;

        const auto g = mjolnir::read_implicit_membrane_potential<real_type>(v);

        BOOST_TEST(g.participants().size() == 3u);
        BOOST_TEST(g.participants().at(0)  == 0u);
        BOOST_TEST(g.participants().at(1)  == 1u);
        BOOST_TEST(g.participants().at(2)  == 5u);
        BOOST_TEST(g.parameters().at(0)    == real_type(  0.0     ), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(1)    == real_type(100.0     ), tolerance<real_type>());
        BOOST_TEST(g.parameters().at(5)    == real_type(100.0     ), tolerance<real_type>());
        BOOST_TEST(g.half_thick()          == real_type(3.14 * 0.5), tolerance<real_type>());
        BOOST_TEST(g.k()                   == real_type(6.28      ), tolerance<real_type>());
        BOOST_TEST(g.bend()                == real_type(9.42      ), tolerance<real_type>());
        BOOST_TEST(g.cutoff()              == real_type(5.0       ), tolerance<real_type>());
    }
}


