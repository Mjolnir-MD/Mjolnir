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
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        output.prefix = "test"
        [[forcefields]]
        # empty
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, toml::table{});
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

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
        const auto v = u8R"(
            [files]
            output.prefix = "test"
            [[forcefields]]
            [[forcefields.local]]
                interaction = "BondLength"
                potential   = "Harmonic"
                topology    = "bond"
                parameters  = []
            [[forcefields.local]]
                interaction = "BondAngle"
                potential   = "Harmonic"
                topology    = "none"
                parameters  = []
            [[forcefields.global]]
                interaction = "Pair"
                potential   = "ExcludedVolume"
                epsilon     = 3.14
                parameters  = [{index = 0, radius = 2.0}]
                ignore.molecule = "Nothing"
                ignore.particles_within = {bond = 3, contact = 1}
                spatial_partition.type = "Naive"
            [[forcefields.global]]
                interaction = "Pair"
                potential   = "LennardJones"
                parameters  = [{index = 0, sigma = 2.0, epsilon = 10.0}]
                ignore.molecule = "Nothing"
                ignore.particles_within = {bond = 3, contact = 1}
                spatial_partition.type = "Naive"
            [[forcefields.external]]
                interaction = "Distance"
                potential   = "ExcludedVolumeWall"
                epsilon     = 3.14
                parameters  = []
                shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
            [[forcefields.external]]
                interaction = "Distance"
                potential   = "LennardJonesWall"
                parameters  = []
                shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, toml::table{});
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(ff.local().size()    == 2u);
        BOOST_TEST(ff.global().size()   == 2u);
        BOOST_TEST(ff.external().size() == 2u);

        // contents are tested in the individual test codes
    }
}
