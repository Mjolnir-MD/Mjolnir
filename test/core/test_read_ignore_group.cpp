#define BOOST_TEST_MODULE "test_read_ignore_group"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_global_potential.hpp>
#include <cmath>

BOOST_AUTO_TEST_CASE(read_ignored_group_default)
{
    mjolnir::LoggerManager::set_default_logger("test_read_ignore_groups");

    const toml::value v = toml::table{}; // empty!
    const auto mol = mjolnir::read_ignored_group(v);

    BOOST_TEST(mol.ignores().empty());

    BOOST_TEST(!mol.is_ignored("protein",  "protein"));
    BOOST_TEST(!mol.is_ignored("protein1", "protein2"));
}

BOOST_AUTO_TEST_CASE(read_ignored_group_intra)
{
    mjolnir::LoggerManager::set_default_logger("test_read_ignore_groups");

    using namespace toml::literals;
    const toml::value v = u8R"(
        group.intra = ["protein"]
    )"_toml;

    const auto mol = mjolnir::read_ignored_group(v);

    BOOST_TEST(mol.ignores().at("protein").size()  == 1u);
    BOOST_TEST(mol.ignores().at("protein").front() == "protein");

    BOOST_TEST( mol.is_ignored("protein",  "protein"));
    BOOST_TEST(!mol.is_ignored("protein1", "protein"));
    BOOST_TEST(!mol.is_ignored("protein",  "protein2"));
}

BOOST_AUTO_TEST_CASE(read_ignored_group_inter)
{
    mjolnir::LoggerManager::set_default_logger("test_read_ignore_groups");

    using namespace toml::literals;
    const toml::value v = u8R"(
        group.inter = [
            ["protein1", "protein2"]
        ]
    )"_toml;

    const auto mol = mjolnir::read_ignored_group(v);

    BOOST_TEST(mol.ignores().at("protein1").size()  == 1u);
    BOOST_TEST(mol.ignores().at("protein1").front() == "protein2");
    BOOST_TEST(mol.ignores().at("protein2").size()  == 1u);
    BOOST_TEST(mol.ignores().at("protein2").front() == "protein1");

    BOOST_TEST(!mol.is_ignored("protein",  "protein"));
    BOOST_TEST( mol.is_ignored("protein1", "protein2"));
}

BOOST_AUTO_TEST_CASE(read_ignored_group_intra_inter)
{
    mjolnir::LoggerManager::set_default_logger("test_read_ignore_groups");

    using namespace toml::literals;
    const toml::value v = u8R"(
        group.intra = ["protein"]
        group.inter = [
            ["protein", "protein1"]
        ]
    )"_toml;

    const auto mol = mjolnir::read_ignored_group(v);

    BOOST_TEST(mol.ignores().at("protein").size() == 2u);
    BOOST_TEST(mol.ignores().at("protein").at(0)  == "protein");
    BOOST_TEST(mol.ignores().at("protein").at(1)  == "protein1");
    BOOST_TEST(mol.ignores().at("protein1").size()  == 1u);
    BOOST_TEST(mol.ignores().at("protein1").front() == "protein");

    BOOST_TEST( mol.is_ignored("protein",  "protein"));
    BOOST_TEST( mol.is_ignored("protein",  "protein1"));
    BOOST_TEST(!mol.is_ignored("protein1", "protein1"));
    BOOST_TEST(!mol.is_ignored("protein1", "protein2"));
}
