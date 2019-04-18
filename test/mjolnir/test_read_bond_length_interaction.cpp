#define BOOST_TEST_MODULE "test_read_bond_length_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_local_interaction.hpp>
#include <test/util/traits.hpp>

BOOST_AUTO_TEST_CASE(read_bond_length_harmonic)
{
    mjolnir::LoggerManager::set_default_logger("test_read_bond_length_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::test::traits<real_type>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "BondLength"
            potential   = "Harmonic"
            topology    = "bond"
            parameters  = [] # empty
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::BondLengthInteraction<
            traits_type, mjolnir::HarmonicPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_bond_length_go_contact)
{
    mjolnir::LoggerManager::set_default_logger("test_read_bond_length_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::test::traits<real_type>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "BondLength"
            potential   = "GoContact"
            topology    = "contact"
            parameters  = [] # empty
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::BondLengthInteraction<
            traits_type, mjolnir::GoContactPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_bond_length_gaussian)
{
    mjolnir::LoggerManager::set_default_logger("test_read_bond_length_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::test::traits<real_type>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "BondLength"
            potential   = "Gaussian"
            topology    = "none"
            parameters  = [] # empty
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::BondLengthInteraction<
            traits_type, mjolnir::GaussianPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}
