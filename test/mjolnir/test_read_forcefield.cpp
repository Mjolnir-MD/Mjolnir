#define BOOST_TEST_MODULE "test_read_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_forcefield.hpp>

BOOST_AUTO_TEST_CASE(read_empty_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        const toml::table v = toml::table{
            {"local",       toml::array{}},
            {"global",      toml::array{}},
            {"external",    toml::array{}}
        };

        const auto ff = mjolnir::read_forcefield_from_table<traits_type>(v, "./");

        BOOST_TEST(ff.local().size()    == 0u);
        BOOST_TEST(ff.global().size()   == 0u);
        BOOST_TEST(ff.external().size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_several_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            [[local]]
                interaction = "BondLength"
                potential   = "Harmonic"
                topology    = "bond"
                parameters  = []
            [[local]]
                interaction = "BondAngle"
                potential   = "Harmonic"
                topology    = "none"
                parameters  = []
            [[global]]
                interaction = "Pair"
                potential   = "ExcludedVolume"
                epsilon     = 3.14
                parameters  = []
                ignore.molecule = "Nothing"
                ignore.particles_within = {bond = 3, contact = 1}
                spatial_partition.type = "Naive"
            [[global]]
                interaction = "Pair"
                potential   = "LennardJones"
                parameters  = []
                ignore.molecule = "Nothing"
                ignore.particles_within = {bond = 3, contact = 1}
                spatial_partition.type = "Naive"
            [[external]]
                interaction = "Distance"
                potential   = "ExcludedVolumeWall"
                epsilon     = 3.14
                parameters  = []
                shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
            [[external]]
                interaction = "Distance"
                potential   = "LennardJonesWall"
                parameters  = []
                shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
        )"_toml);

        const auto ff = mjolnir::read_forcefield_from_table<traits_type>(v, "./");
        BOOST_TEST(ff.local().size()    == 2u);
        BOOST_TEST(ff.global().size()   == 2u);
        BOOST_TEST(ff.external().size() == 2u);

        // contents are tested in the individual test codes
    }
}
