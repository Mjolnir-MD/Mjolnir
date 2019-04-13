#define BOOST_TEST_MODULE "test_read_harmonic_restraint_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_external_potential.hpp>

BOOST_AUTO_TEST_CASE(read_harmonic_restraint_double)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic_restraint.log");

    using real_type = double;
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
        BOOST_TEST(h.parameters().at(0).first  == 10.0);
        BOOST_TEST(h.parameters().at(0).second ==  0.0);
        BOOST_TEST(h.parameters().at(1).first  == 20.0);
        BOOST_TEST(h.parameters().at(1).second ==  1.0);
    }
}

BOOST_AUTO_TEST_CASE(read_harmonic_restraint_float)
{
    mjolnir::LoggerManager::set_default_logger("test_read_harmonic_restraint.log");
    using real_type = float;
    constexpr real_type tol = 1e-4; // for conversion between double <-> float

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
        BOOST_TEST(h.parameters().at(0).first  == 10.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(h.parameters().at(0).second ==  0.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(h.parameters().at(1).first  == 20.0f, boost::test_tools::tolerance(tol));
        BOOST_TEST(h.parameters().at(1).second ==  1.0f, boost::test_tools::tolerance(tol));
    }
}
