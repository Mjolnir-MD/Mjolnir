#define BOOST_TEST_MODULE "test_read_hybrid_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_forcefield.hpp>

BOOST_AUTO_TEST_CASE(read_hybrid_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_hybrid_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        output.prefix = "test"
        [simulator]
        forcefields.type = "Hybrid"
        forcefields.lambda = 0.5
        [[forcefields]]
            name = "open"
          [[forcefields.local]]
            interaction = "BondLength"
            potential   = "Harmonic"
            topology    = "bond"
            parameters  = []
          [[forcefields.global]]
            interaction = "Pair"
            potential   = "ExcludedVolume"
            epsilon     = 3.14
            parameters  = []
            ignore.molecule = "Nothing"
            ignore.particles_within = {bond = 3, contact = 1}
            spatial_partition.type = "Naive"
          [[forcefields.external]]
            interaction = "Distance"
            potential   = "ExcludedVolumeWall"
            epsilon     = 3.14
            parameters  = []
            shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
        [[forcefields]]
            name = "close"
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
            parameters  = []
            ignore.molecule = "Nothing"
            ignore.particles_within = {bond = 3, contact = 1}
            spatial_partition.type = "Naive"
          [[forcefields.global]]
            interaction = "Pair"
            potential   = "LennardJones"
            parameters  = []
            ignore.molecule = "Nothing"
            ignore.particles_within = {bond = 3, contact = 1}
            spatial_partition.type = "Naive"
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, v.at("simulator"));
        const auto ffp = dynamic_cast<mjolnir::HybridForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto& ff = *ffp;

        BOOST_TEST(ff.lambda() == 0.5);

        const auto ff1 = dynamic_cast<mjolnir::ForceField<traits_type> const*>(ff.ff1().get());
        const auto ff2 = dynamic_cast<mjolnir::ForceField<traits_type> const*>(ff.ff2().get());
        BOOST_REQUIRE(static_cast<bool>(ff1));
        BOOST_REQUIRE(static_cast<bool>(ff2));

        BOOST_TEST(ff1->local().size()    == 1u);
        BOOST_TEST(ff1->global().size()   == 1u);
        BOOST_TEST(ff1->external().size() == 1u);

        BOOST_TEST(ff2->local().size()    == 2u);
        BOOST_TEST(ff2->global().size()   == 2u);
        BOOST_TEST(ff2->external().size() == 0u);
    }
}
