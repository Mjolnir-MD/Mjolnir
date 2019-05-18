#define BOOST_TEST_MODULE "test_read_external_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_external_forcefield.hpp>

#include <typeindex>
#include <typeinfo>

BOOST_AUTO_TEST_CASE(read_empty_external_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_external_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        const toml::array v{};
        const auto ff = mjolnir::read_external_forcefield<traits_type>(v, "./");
        BOOST_TEST(ff.empty());
        BOOST_TEST(ff.size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_external_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_external_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::array v = {u8R"(
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
        )"_toml};

        const auto ff = mjolnir::read_external_forcefield<traits_type>(v, "./");
        BOOST_TEST(!ff.empty());
        BOOST_TEST(ff.size() == 1u);

        const auto& interaction_ptr = *ff.begin();
        BOOST_TEST(static_cast<bool>(interaction_ptr));

        const auto derived_ptr  = dynamic_cast<mjolnir::ExternalDistanceInteraction<
            traits_type, mjolnir::ExcludedVolumeWallPotential<real_type>,
            mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveXDirection>
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
        const toml::array v = {u8R"(
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
        )"_toml, u8R"(
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
        )"_toml};

        const auto ff = mjolnir::read_external_forcefield<traits_type>(v, "./");
        BOOST_TEST(!ff.empty());
        BOOST_TEST(ff.size() == 2u);

        using exv_interaction = mjolnir::ExternalDistanceInteraction<
            traits_type, mjolnir::ExcludedVolumeWallPotential<real_type>,
            mjolnir::AxisAlignedPlane<traits_type, mjolnir::NegativeXDirection>
            >;
        using lj_interaction = mjolnir::ExternalDistanceInteraction<
            traits_type, mjolnir::LennardJonesWallPotential<real_type>,
            mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveXDirection>
            >;

        std::map<std::type_index, bool> found;
        found[typeid(exv_interaction)] = false;
        found[typeid( lj_interaction)] = false;

        for(const auto& interaction_ptr : ff)
        {
            BOOST_TEST(static_cast<bool>(interaction_ptr));
            const auto& deref = *interaction_ptr;
            found[typeid(deref)] = true;
        }
        BOOST_TEST(found[typeid(exv_interaction)]);
        BOOST_TEST(found[typeid( lj_interaction)]);
    }
}
