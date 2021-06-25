#define BOOST_TEST_MODULE "test_read_global_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_forcefield.hpp>

#include <typeindex>
#include <typeinfo>

BOOST_AUTO_TEST_CASE(read_empty_global_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            [files]
            output.prefix = "test"
            [[forcefields]]
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, toml::table{});
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(ff.global().empty());
        BOOST_TEST(ff.global().size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_global_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            [files]
            output.prefix = "test"
            [[forcefields]]
            [[forcefields.global]]
            interaction                     = "Pair"
            potential                       = "ExcludedVolume"
            spatial_partition.type          = "Naive"
            epsilon                         = 3.14
            ignore.molecule                 = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = []
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, toml::table{});
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(!ff.global().empty());
        BOOST_TEST( ff.global().size() == 1u);

        const auto& interaction_ptr = *(ff.global().begin());
        BOOST_TEST(static_cast<bool>(interaction_ptr));

        using exv_interaction = mjolnir::GlobalPairInteraction<
                traits_type, mjolnir::ExcludedVolumePotential<real_type>
            >;

        const auto derived_ptr  = dynamic_cast<exv_interaction*>(interaction_ptr.get());
        BOOST_TEST(static_cast<bool>(derived_ptr));
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
            [[forcefields.global]]
            interaction                     = "Pair"
            potential                       = "ExcludedVolume"
            spatial_partition.type          = "Naive"
            epsilon                         = 3.14
            ignore.molecule                 = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = []
            [[forcefields.global]]
            interaction                     = "Pair"
            potential                       = "LennardJones"
            spatial_partition.type          = "Naive"
            ignore.molecule                 = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = []
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, toml::table{});
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(!ff.global().empty());
        BOOST_TEST( ff.global().size() == 2u);

        using exv_interaction = mjolnir::GlobalPairInteraction<
                traits_type, mjolnir::ExcludedVolumePotential<real_type>
            >;
        using lj_interaction  = mjolnir::GlobalPairInteraction<
                traits_type, mjolnir::LennardJonesPotential<real_type>
            >;

        std::map<std::type_index, bool> found;
        found[typeid(exv_interaction)] = false;
        found[typeid(lj_interaction)]  = false;

        for(const auto& interaction_ptr : ff.global())
        {
            BOOST_TEST(static_cast<bool>(interaction_ptr));
            const auto& deref = *interaction_ptr;
            found[typeid(deref)] = true;
        }
        BOOST_TEST(found[typeid(exv_interaction)]);
        BOOST_TEST(found[typeid(lj_interaction)]);
    }
}
