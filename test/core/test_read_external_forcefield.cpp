#define BOOST_TEST_MODULE "test_read_external_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_forcefield.hpp>

#include <typeindex>
#include <typeinfo>

BOOST_AUTO_TEST_CASE(read_empty_external_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_external_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            [files]
            output.prefix = "test"
            [[forcefields]]
        )"_toml;
        const auto ffb = mjolnir::read_forcefield<traits_type>(v, 0);
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(ff.external().empty());
        BOOST_TEST(ff.external().size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_external_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_external_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            [files]
            output.prefix = "test"
            [[forcefields]]
            [[forcefields.external]]
            interaction    = "Distance"
            potential      = "ExcludedVolumeWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "+X"
            shape.position = 1.0
            shape.margin   = 0.5
            epsilon        = 3.14
            parameters     = [
                {index = 0, radius = 2.0},
                {index = 1, radius = 2.0},
            ]
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, 0);
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(!ff.external().empty());
        BOOST_TEST(ff.external().size() == 1u);

        const auto& interaction_ptr = *(ff.external().begin());
        BOOST_TEST(static_cast<bool>(interaction_ptr));

        const auto derived_ptr  = dynamic_cast<mjolnir::ExternalDistanceInteraction<
            traits_type, mjolnir::ExcludedVolumeWallPotential<real_type>,
            mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveXDirection<traits_type>>
            >*>(interaction_ptr.get());
        BOOST_TEST(static_cast<bool>(derived_ptr));
    }
}

BOOST_AUTO_TEST_CASE(read_several_external_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_external_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            [files]
            output.prefix = "test"
            [[forcefields]]
            [[forcefields.external]]
            interaction    = "Distance"
            potential      = "ExcludedVolumeWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "-X"
            shape.position = 1.0
            shape.margin   = 0.5
            epsilon        = 3.14
            parameters     = [
                {index = 0, radius = 2.0},
                {index = 1, radius = 2.0},
            ]
            [[forcefields.external]]
            interaction    = "Distance"
            potential      = "LennardJonesWall"
            shape.name     = "AxisAlignedPlane"
            shape.axis     = "+X"
            shape.position = 1.0
            shape.margin   = 0.5
            parameters     = [
                {index = 0, sigma = 2.0, epsilon = 2.0},
                {index = 1, sigma = 2.0, epsilon = 2.0},
            ]
        )"_toml;

        const auto ffb = mjolnir::read_forcefield<traits_type>(v, 0);
        const auto ffp = dynamic_cast<mjolnir::ForceField<traits_type>*>(ffb.get());
        BOOST_REQUIRE(ffp);
        const auto ff = *ffp;

        BOOST_TEST(!ff.external().empty());
        BOOST_TEST(ff.external().size() == 2u);

        using exv_interaction = mjolnir::ExternalDistanceInteraction<
            traits_type, mjolnir::ExcludedVolumeWallPotential<real_type>,
            mjolnir::AxisAlignedPlane<traits_type, mjolnir::NegativeXDirection<traits_type>>
            >;
        using lj_interaction = mjolnir::ExternalDistanceInteraction<
            traits_type, mjolnir::LennardJonesWallPotential<real_type>,
            mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveXDirection<traits_type>>
            >;

        std::map<std::type_index, bool> found;
        found[typeid(exv_interaction)] = false;
        found[typeid( lj_interaction)] = false;

        for(const auto& interaction_ptr : ff.external())
        {
            BOOST_TEST(static_cast<bool>(interaction_ptr));
            const auto& deref = *interaction_ptr;
            found[typeid(deref)] = true;
        }
        BOOST_TEST(found[typeid(exv_interaction)]);
        BOOST_TEST(found[typeid( lj_interaction)]);
    }
}
