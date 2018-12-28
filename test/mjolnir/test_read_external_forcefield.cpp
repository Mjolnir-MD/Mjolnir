#define BOOST_TEST_MODULE "test_read_external_forcefield"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_external_forcefield.hpp>

#include <typeindex>
#include <typeinfo>

BOOST_AUTO_TEST_CASE(read_empty_external_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_external_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const toml::array v{};
        const auto ff = mjolnir::read_external_forcefield<traits_type>(v, "./");
        BOOST_TEST(ff.empty());
        BOOST_TEST(ff.size() == 0);
    }
}

BOOST_AUTO_TEST_CASE(read_external_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_external_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const toml::array v{toml::table{
            {"interaction",       toml::value("Distance")},
            {"potential",         toml::value("ExcludedVolumeWall")},
            {"shape", toml::value(toml::table{
                    {"name",     toml::value("AxisAlignedPlane")},
                    {"axis",     toml::value("+X")},
                    {"position", toml::value(1.0)},
                    {"margin",   toml::value(0.5)}
            })},
            {"epsilon",           toml::value(3.14)},
            {"parameters",        toml::value(toml::array{
            })}
        }};

        const auto ff = mjolnir::read_external_forcefield<traits_type>(v, "./");
        BOOST_TEST(!ff.empty());
        BOOST_TEST(ff.size() == 1);

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
    constexpr real_type tol = 1e-8;
    {
        const toml::array v{toml::table{
                {"interaction",       toml::value("Distance")},
                {"potential",         toml::value("ExcludedVolumeWall")},
                {"shape", toml::value(toml::table{
                        {"name",     toml::value("AxisAlignedPlane")},
                        {"axis",     toml::value("+X")},
                        {"position", toml::value(1.0)},
                        {"margin",   toml::value(0.5)}
                })},
                {"epsilon",           toml::value(3.14)},
                {"parameters",        toml::value(toml::array{})}
            }, toml::table{
                {"interaction",       toml::value("Distance")},
                {"potential",         toml::value("LennardJonesWall")},
                {"shape", toml::value(toml::table{
                        {"name",     toml::value("AxisAlignedPlane")},
                        {"axis",     toml::value("+X")},
                        {"position", toml::value(1.0)},
                        {"margin",   toml::value(0.5)}
                })},
                {"parameters",        toml::value(toml::array{})}
            }
        };

        const auto ff = mjolnir::read_external_forcefield<traits_type>(v, "./");
        BOOST_TEST(!ff.empty());
        BOOST_TEST(ff.size() == 2);

        using exv_interaction = mjolnir::ExternalDistanceInteraction<
            traits_type, mjolnir::ExcludedVolumeWallPotential<real_type>,
            mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveXDirection>
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
            found[typeid(*interaction_ptr)] = true;
        }
        BOOST_TEST(found[typeid(exv_interaction)]);
        BOOST_TEST(found[typeid( lj_interaction)]);
    }
}
