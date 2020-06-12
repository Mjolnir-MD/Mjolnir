#define BOOST_TEST_MODULE "test_read_local_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_forcefield.hpp>

#include <typeindex>
#include <typeinfo>

BOOST_AUTO_TEST_CASE(read_empty_local_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_local_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        output.prefix = "name"
        [[forcefields]]
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, toml::value(toml::table{}));
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(ff.local().empty());
        BOOST_TEST(ff.local().size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_local_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_local_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        # empty
        [[forcefields]]
        [[forcefields.local]]
        interaction = "BondAngle"
        potential   = "Harmonic"
        topology    = "none"
        parameters  = []
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, toml::value(toml::table{}));
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(!ff.local().empty());
        BOOST_TEST(ff.local().size() == 1u);

        const auto& interaction_ptr = *(ff.local().begin());
        BOOST_TEST(static_cast<bool>(interaction_ptr));

        const auto bond_angle_ptr  = dynamic_cast<mjolnir::BondAngleInteraction<
            traits_type, mjolnir::HarmonicPotential<real_type>>*
            >(interaction_ptr.get());
        BOOST_TEST(static_cast<bool>(bond_angle_ptr));
    }
}

BOOST_AUTO_TEST_CASE(read_several_local_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_local_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
        [files]
        # empty
        [[forcefields]]
        [[forcefields.local]]
        interaction = "BondAngle"
        potential   = "Harmonic"
        topology    = "none"
        parameters  = []
        [[forcefields.local]]
        interaction = "BondLength"
        potential   = "Harmonic"
        topology    = "bond"
        parameters  = []
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, toml::value(toml::table{}));
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(!ff.local().empty());
        BOOST_TEST( ff.local().size() == 2u);

        using bond_length_interaction = mjolnir::BondLengthInteraction<
            traits_type, mjolnir::HarmonicPotential<real_type>>;
        using bond_angle_interaction  = mjolnir::BondAngleInteraction<
            traits_type, mjolnir::HarmonicPotential<real_type>>;

        std::map<std::type_index, bool> found;
        found[typeid(bond_length_interaction)] = false;
        found[typeid(bond_angle_interaction)]  = false;

        for(const auto& interaction_ptr : ff.local())
        {
            BOOST_TEST(static_cast<bool>(interaction_ptr));
            const auto& deref = *interaction_ptr;
            found[typeid(deref)] = true;
        }
        BOOST_TEST(found[typeid(bond_length_interaction)]);
        BOOST_TEST(found[typeid(bond_angle_interaction)]);
    }
}
