#define BOOST_TEST_MODULE "test_read_multiple_basin_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/input/read_forcefield.hpp>

BOOST_AUTO_TEST_CASE(read_empty_forcefield_2basin)
{
    mjolnir::LoggerManager::set_default_logger("test_read_multiple_basin_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        output.prefix = "test"
        [simulator]
        forcefields.type = "MultipleBasin"
        forcefields.units = [
            {basins=["open", "close"], dVs=[0.0, 10.0], delta = 60.0}
        ]
        [[forcefields]]
        name = "open"
        # empty
        [[forcefields]]
        name = "close"
        # empty
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, v.at("simulator"));
        const auto ffp = dynamic_cast<mjolnir::MultipleBasinForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto& ff = *ffp;

        BOOST_TEST(ff.common_local().size()    == 0u);
        BOOST_TEST(ff.common_global().size()   == 0u);
        BOOST_TEST(ff.common_external().size() == 0u);

        const auto& units = ff.units();
        BOOST_REQUIRE(units.size() == 1u);

        const auto  unit1_ptr = dynamic_cast<mjolnir::MultipleBasin2BasinUnit<traits_type>*>(units.at(0).get());
        BOOST_REQUIRE(unit1_ptr);
        const auto& unit1 = *unit1_ptr;

        BOOST_TEST(unit1.delta() == -60.0);
        BOOST_TEST(unit1.dV1()   ==   0.0);
        BOOST_TEST(unit1.dV2()   ==  10.0);

        BOOST_TEST(unit1.name1() == "open");
        BOOST_TEST(unit1.name2() == "close");

        BOOST_TEST(unit1.local1().size()    == 0u);
        BOOST_TEST(unit1.local2().size()    == 0u);
        BOOST_TEST(unit1.global1().size()   == 0u);
        BOOST_TEST(unit1.global2().size()   == 0u);
        BOOST_TEST(unit1.external1().size() == 0u);
        BOOST_TEST(unit1.external2().size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_empty_forcefield_3basin)
{
    mjolnir::LoggerManager::set_default_logger("test_read_multiple_basin_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        output.prefix = "test"
        [simulator]
        forcefields.type = "MultipleBasin"
        forcefields.units = [
            {basins=["open", "mid", "close"], dVs=[0.0, -5.0, 10.0], delta.open-mid = 60.0, delta.mid-close=50.0, delta.close-open = 100.0}
        ]
        [[forcefields]]
        name = "open"
        # empty
        [[forcefields]]
        name = "mid"
        # empty
        [[forcefields]]
        name = "close"
        # empty
        [[forcefields]]
        name = "common"
        # empty
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, v.at("simulator"));
        const auto ffp = dynamic_cast<mjolnir::MultipleBasinForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto& ff = *ffp;

        BOOST_TEST(ff.common_local().size()    == 0u);
        BOOST_TEST(ff.common_global().size()   == 0u);
        BOOST_TEST(ff.common_external().size() == 0u);

        const auto& units = ff.units();
        BOOST_REQUIRE(units.size() == 1u);

        const auto  unit1_ptr = dynamic_cast<mjolnir::MultipleBasin3BasinUnit<traits_type>*>(units.at(0).get());
        BOOST_REQUIRE(unit1_ptr);
        const auto& unit1 = *unit1_ptr;

        BOOST_TEST(unit1.delta12() ==  -60.0);
        BOOST_TEST(unit1.delta23() ==  -50.0);
        BOOST_TEST(unit1.delta31() == -100.0);
        BOOST_TEST(unit1.dV1()     ==   0.0);
        BOOST_TEST(unit1.dV2()     ==  -5.0);
        BOOST_TEST(unit1.dV3()     ==  10.0);

        BOOST_TEST(unit1.name1() == "open");
        BOOST_TEST(unit1.name2() == "mid");
        BOOST_TEST(unit1.name3() == "close");

        BOOST_TEST(unit1.local1().size()    == 0u);
        BOOST_TEST(unit1.local2().size()    == 0u);
        BOOST_TEST(unit1.local3().size()    == 0u);
        BOOST_TEST(unit1.global1().size()   == 0u);
        BOOST_TEST(unit1.global2().size()   == 0u);
        BOOST_TEST(unit1.global3().size()   == 0u);
        BOOST_TEST(unit1.external1().size() == 0u);
        BOOST_TEST(unit1.external2().size() == 0u);
        BOOST_TEST(unit1.external3().size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_several_forcefield_2basin)
{
    mjolnir::LoggerManager::set_default_logger("test_read_multiple_basin_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        output.prefix = "test"
        [simulator]
        forcefields.type = "MultipleBasin"
        forcefields.units = [
            {basins=["open", "close"], dVs=[0.0, 10.0], delta = 60.0}
        ]
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
        [[forcefields]]
            name = "common"
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
          [[forcefields.local]]
            interaction = "DihedralAngle"
            potential   = "Gaussian"
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
          [[forcefields.global]]
            interaction = "Pair"
            potential   = "DebyeHuckel"
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
          [[forcefields.external]]
            interaction = "Distance"
            potential   = "ExcludedVolumeWall"
            epsilon     = 3.14
            parameters  = []
            shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}

        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, v.at("simulator"));
        const auto ffp = dynamic_cast<mjolnir::MultipleBasinForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto& ff = *ffp;

        BOOST_TEST(ff.common_local().size()    == 3u);
        BOOST_TEST(ff.common_global().size()   == 3u);
        BOOST_TEST(ff.common_external().size() == 2u);

        const auto& units = ff.units();
        BOOST_REQUIRE(units.size() == 1u);

        const auto  unit1_ptr = dynamic_cast<mjolnir::MultipleBasin2BasinUnit<traits_type>*>(units.at(0).get());
        BOOST_REQUIRE(unit1_ptr);
        const auto& unit1 = *unit1_ptr;

        BOOST_TEST(unit1.delta() == -60.0);
        BOOST_TEST(unit1.dV1()   ==   0.0);
        BOOST_TEST(unit1.dV2()   ==  10.0);

        BOOST_TEST(unit1.name1() == "open");
        BOOST_TEST(unit1.name2() == "close");

        BOOST_TEST(unit1.local1().size()    == 1u);
        BOOST_TEST(unit1.global1().size()   == 1u);
        BOOST_TEST(unit1.external1().size() == 1u);

        BOOST_TEST(unit1.local2().size()    == 2u);
        BOOST_TEST(unit1.global2().size()   == 2u);
        BOOST_TEST(unit1.external2().size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_several_forcefield_3basin)
{
    mjolnir::LoggerManager::set_default_logger("test_read_multiple_basin_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        output.prefix = "test"
        [simulator]
        [simulator.forcefields]
        type = "MultipleBasin"
        [[simulator.forcefields.units]]
        basins = ["open", "mid", "close"]
        dVs    = [   0.0,  -5.0,    10.0]
        delta.open-mid   = 60.0
        delta.mid-close  = 50.0
        delta.close-open = 100.0

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
            name = "mid"
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
          [[forcefields.external]]
            interaction = "Distance"
            potential   = "ExcludedVolumeWall"
            epsilon     = 3.14
            parameters  = []
            shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
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
          [[forcefields.external]]
            interaction = "Distance"
            potential   = "ExcludedVolumeWall"
            epsilon     = 3.14
            parameters  = []
            shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
          [[forcefields.external]]
            interaction = "Distance"
            potential   = "ExcludedVolumeWall"
            epsilon     = 3.14
            parameters  = []
            shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}

        [[forcefields]]
            name = "common"
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
          [[forcefields.local]]
            interaction = "DihedralAngle"
            potential   = "Gaussian"
            topology    = "none"
            parameters  = []
          [[forcefields.local]]
            interaction = "DihedralAngle"
            potential   = "Gaussian"
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
          [[forcefields.global]]
            interaction = "Pair"
            potential   = "DebyeHuckel"
            parameters  = []
            ignore.molecule = "Nothing"
            ignore.particles_within = {bond = 3, contact = 1}
            spatial_partition.type = "Naive"
          [[forcefields.global]]
            interaction = "Pair"
            potential   = "DebyeHuckel"
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
          [[forcefields.external]]
            interaction = "Distance"
            potential   = "ExcludedVolumeWall"
            epsilon     = 3.14
            parameters  = []
            shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
          [[forcefields.external]]
            interaction = "Distance"
            potential   = "ExcludedVolumeWall"
            epsilon     = 3.14
            parameters  = []
            shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}
          [[forcefields.external]]
            interaction = "Distance"
            potential   = "ExcludedVolumeWall"
            epsilon     = 3.14
            parameters  = []
            shape = {name = "AxisAlignedPlane", axis="+X", position=1.0, margin=0.5}

        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, v.at("simulator"));
        const auto ffp = dynamic_cast<mjolnir::MultipleBasinForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto& ff = *ffp;

        BOOST_TEST(ff.common_local().size()    == 4u);
        BOOST_TEST(ff.common_global().size()   == 4u);
        BOOST_TEST(ff.common_external().size() == 4u);

        const auto& units = ff.units();
        BOOST_REQUIRE(units.size() == 1u);

        const auto  unit1_ptr = dynamic_cast<mjolnir::MultipleBasin3BasinUnit<traits_type>*>(units.at(0).get());
        BOOST_REQUIRE(unit1_ptr);
        const auto& unit1 = *unit1_ptr;

        BOOST_TEST(unit1.delta12() ==  -60.0);
        BOOST_TEST(unit1.delta23() ==  -50.0);
        BOOST_TEST(unit1.delta31() == -100.0);
        BOOST_TEST(unit1.dV1()     ==   0.0);
        BOOST_TEST(unit1.dV2()     ==  -5.0);
        BOOST_TEST(unit1.dV3()     ==  10.0);

        BOOST_TEST(unit1.name1() == "open");
        BOOST_TEST(unit1.name2() == "mid");
        BOOST_TEST(unit1.name3() == "close");

        BOOST_TEST(unit1.local1().size()    == 1u);
        BOOST_TEST(unit1.global1().size()   == 1u);
        BOOST_TEST(unit1.external1().size() == 1u);

        BOOST_TEST(unit1.local2().size()    == 2u);
        BOOST_TEST(unit1.global2().size()   == 2u);
        BOOST_TEST(unit1.external2().size() == 2u);

        BOOST_TEST(unit1.local3().size()    == 3u);
        BOOST_TEST(unit1.global3().size()   == 3u);
        BOOST_TEST(unit1.external3().size() == 3u);
    }
}

BOOST_AUTO_TEST_CASE(read_empty_forcefield_2basins_2units)
{
    mjolnir::LoggerManager::set_default_logger("test_read_multiple_basin_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        output.prefix = "test"
        [simulator]
        forcefields.type = "MultipleBasin"
        forcefields.units = [
            {basins=["domain1-open", "domain1-close"], dVs=[0.0, 10.0], delta = 60.0},
            {basins=["domain2-open", "domain2-close"], dVs=[0.0, 15.0], delta = 50.0}
        ]
        [[forcefields]]
        name = "domain1-open"
        # empty
        [[forcefields]]
        name = "domain1-close"
        # empty
        [[forcefields]]
        name = "domain2-open"
        # empty
        [[forcefields]]
        name = "domain2-close"
        # empty
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, v.at("simulator"));
        const auto ffp = dynamic_cast<mjolnir::MultipleBasinForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto& ff = *ffp;

        BOOST_TEST(ff.common_local().size()    == 0u);
        BOOST_TEST(ff.common_global().size()   == 0u);
        BOOST_TEST(ff.common_external().size() == 0u);

        const auto& units = ff.units();
        BOOST_REQUIRE(units.size() == 2u);

        const auto  unit1_ptr = dynamic_cast<mjolnir::MultipleBasin2BasinUnit<traits_type>*>(units.at(0).get());
        BOOST_REQUIRE(unit1_ptr);
        const auto& unit1 = *unit1_ptr;

        BOOST_TEST(unit1.delta() == -60.0);
        BOOST_TEST(unit1.dV1()   ==   0.0);
        BOOST_TEST(unit1.dV2()   ==  10.0);

        BOOST_TEST(unit1.name1() == "domain1-open");
        BOOST_TEST(unit1.name2() == "domain1-close");

        BOOST_TEST(unit1.local1().size()    == 0u);
        BOOST_TEST(unit1.local2().size()    == 0u);
        BOOST_TEST(unit1.global1().size()   == 0u);
        BOOST_TEST(unit1.global2().size()   == 0u);
        BOOST_TEST(unit1.external1().size() == 0u);
        BOOST_TEST(unit1.external2().size() == 0u);

        const auto  unit2_ptr = dynamic_cast<mjolnir::MultipleBasin2BasinUnit<traits_type>*>(units.at(1).get());
        BOOST_REQUIRE(unit2_ptr);
        const auto& unit2 = *unit2_ptr;

        BOOST_TEST(unit2.delta() == -50.0);
        BOOST_TEST(unit2.dV1()   ==   0.0);
        BOOST_TEST(unit2.dV2()   ==  15.0);

        BOOST_TEST(unit2.name1() == "domain2-open");
        BOOST_TEST(unit2.name2() == "domain2-close");

        BOOST_TEST(unit2.local1().size()    == 0u);
        BOOST_TEST(unit2.local2().size()    == 0u);
        BOOST_TEST(unit2.global1().size()   == 0u);
        BOOST_TEST(unit2.global2().size()   == 0u);
        BOOST_TEST(unit2.external1().size() == 0u);
        BOOST_TEST(unit2.external2().size() == 0u);

    }
}
