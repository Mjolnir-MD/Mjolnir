#define BOOST_TEST_MODULE "test_read_local_forcefield"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_local_forcefield.hpp>

#include <typeindex>
#include <typeinfo>

BOOST_AUTO_TEST_CASE(read_empty_local_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_local_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const std::vector<toml::table> v{};
        const auto ff = mjolnir::read_local_forcefield<traits_type>(v, "./");
        BOOST_TEST(ff.empty());
        BOOST_TEST(ff.size() == 0);
    }
}

BOOST_AUTO_TEST_CASE(read_local_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_local_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const std::vector<toml::table> v{toml::table{
                {"interaction", toml::value("BondAngle")},
                {"potential",   toml::value("Harmonic")},
                {"topology",    toml::value("none")},
                {"parameters",  toml::value(toml::array(/*empty*/))}
            }
        };

        const auto lff = mjolnir::read_local_forcefield<traits_type>(v, "./");
        BOOST_TEST(!lff.empty());
        BOOST_TEST(lff.size() == 1);

        const auto& interaction_ptr = *lff.begin();
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
    constexpr real_type tol = 1e-8;
    {
        const std::vector<toml::table> v{toml::table{
                {"interaction", toml::value("BondAngle")},
                {"potential",   toml::value("Harmonic")},
                {"topology",    toml::value("none")},
                {"parameters",  toml::value(toml::array(/*empty*/))}
            }, toml::table{
                {"interaction", toml::value("BondLength")},
                {"potential",   toml::value("Harmonic")},
                {"topology",    toml::value("bond")},
                {"parameters",  toml::value(toml::array(/*empty*/))}
            }
        };

        const auto lff = mjolnir::read_local_forcefield<traits_type>(v, "./");
        BOOST_TEST(!lff.empty());
        BOOST_TEST(lff.size() == 2);

        using bond_length_interaction = mjolnir::BondLengthInteraction<
            traits_type, mjolnir::HarmonicPotential<real_type>>;
        using bond_angle_interaction  = mjolnir::BondAngleInteraction<
            traits_type, mjolnir::HarmonicPotential<real_type>>;

        std::map<std::type_index, bool> found;
        found[typeid(bond_length_interaction)] = false;
        found[typeid(bond_angle_interaction)]  = false;

        for(const auto& interaction_ptr : lff)
        {
            BOOST_TEST(static_cast<bool>(interaction_ptr));
            found[typeid(*interaction_ptr)] = true;
        }
        BOOST_TEST(found[typeid(bond_length_interaction)]);
        BOOST_TEST(found[typeid(bond_angle_interaction)]);
    }
}
