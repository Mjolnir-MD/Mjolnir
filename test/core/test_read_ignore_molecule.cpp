#define BOOST_TEST_MODULE "test_read_ignore_molecule"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_global_potential.hpp>
#include <cmath>

BOOST_AUTO_TEST_CASE(read_ignored_molecule_default)
{
    mjolnir::LoggerManager::set_default_logger("test_read_ignore_molecules");

    toml::value v = toml::table{}; // empty!
    const auto mol = mjolnir::read_ignored_molecule(v);

    BOOST_TEST(mol.name() == "Nothing");
}

BOOST_AUTO_TEST_CASE(read_ignored_molecule_nothing)
{
    mjolnir::LoggerManager::set_default_logger("test_read_ignore_molecules");

    using namespace toml::literals;
    const toml::value v = u8R"(
        ignore.molecule = "Nothing"
    )"_toml;

    const auto mol = mjolnir::read_ignored_molecule(v);

    BOOST_TEST(mol.name() == "Nothing");
}

BOOST_AUTO_TEST_CASE(read_ignored_molecule_self)
{
    mjolnir::LoggerManager::set_default_logger("test_read_ignore_molecules");

    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            ignore.molecule = "Self"
        )"_toml;

        const auto mol = mjolnir::read_ignored_molecule(v);

        BOOST_TEST(mol.name() == "Self");
    }
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            ignore.molecule = "Intra"
        )"_toml;

        const auto mol = mjolnir::read_ignored_molecule(v);

        BOOST_TEST(mol.name() == "Self");
    }
}

BOOST_AUTO_TEST_CASE(read_ignored_molecule_others)
{
    mjolnir::LoggerManager::set_default_logger("test_read_ignore_molecules");

    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            ignore.molecule = "Others"
        )"_toml;
        const auto mol = mjolnir::read_ignored_molecule(v);

        BOOST_TEST(mol.name() == "Others");
    }
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            ignore.molecule = "Inter"
        )"_toml;
        const auto mol = mjolnir::read_ignored_molecule(v);

        BOOST_TEST(mol.name() == "Others");
    }
}
