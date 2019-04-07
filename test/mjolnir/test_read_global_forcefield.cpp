#define BOOST_TEST_MODULE "test_read_global_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_global_forcefield.hpp>

#include <typeindex>
#include <typeinfo>

BOOST_AUTO_TEST_CASE(read_empty_global_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        const toml::array v{};
        const auto ff = mjolnir::read_global_forcefield<traits_type>(v, "./");
        BOOST_TEST(ff.empty());
        BOOST_TEST(ff.size() == 0u);
    }
}

BOOST_AUTO_TEST_CASE(read_global_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        const toml::array v{toml::table{
            {"interaction",       toml::value("Pair")},
            {"potential",         toml::value("ExcludedVolume")},
            {"spatial_partition", toml::value(toml::table{
                        {"type", toml::value("Naive")}
            })},
            {"epsilon",           toml::value(3.14)},
            {"ignore",            toml::value(toml::table{
                {"molecule",         toml::value("Nothing")},
                {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
            })},
            {"parameters",        toml::value(toml::array())}
            }
        };

        const auto ff = mjolnir::read_global_forcefield<traits_type>(v, "./");
        BOOST_TEST(!ff.empty());
        BOOST_TEST(ff.size() == 1u);

        const auto& interaction_ptr = *ff.begin();
        BOOST_TEST(static_cast<bool>(interaction_ptr));

        using exv_interaction = mjolnir::GlobalPairInteraction<
                traits_type, mjolnir::ExcludedVolumePotential<real_type>,
                mjolnir::NaivePairCalculation<traits_type,
                    typename mjolnir::ExcludedVolumePotential<real_type>::parameter_type>
            >;

        const auto derived_ptr  = dynamic_cast<exv_interaction*>(interaction_ptr.get());
        BOOST_TEST(static_cast<bool>(derived_ptr));
    }
}

BOOST_AUTO_TEST_CASE(read_several_global_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    {
        const toml::array v{toml::table{
                {"interaction",       toml::value("Pair")},
                {"potential",         toml::value("ExcludedVolume")},
                {"spatial_partition", toml::value(toml::table{
                            {"type", toml::value("Naive")}
                })},
                {"epsilon",           toml::value(3.14)},
                {"ignore",            toml::value(toml::table{
                    {"molecule",         toml::value("Nothing")},
                    {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
                })},
                {"parameters",        toml::value(toml::array())}
            }, toml::table{
                {"interaction",       toml::value("Pair")},
                {"potential",         toml::value("LennardJones")},
                {"spatial_partition", toml::value(toml::table{
                            {"type", toml::value("Naive")}
                })},
                {"ignore",            toml::value(toml::table{
                    {"molecule",         toml::value("Nothing")},
                    {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
                })},
                {"parameters",        toml::value(toml::array{})}
            }
        };

        const auto ff = mjolnir::read_global_forcefield<traits_type>(v, "./");
        BOOST_TEST(!ff.empty());
        BOOST_TEST(ff.size() == 2u);

        using exv_interaction = mjolnir::GlobalPairInteraction<
                traits_type, mjolnir::ExcludedVolumePotential<real_type>,
                mjolnir::NaivePairCalculation<traits_type,
                    typename mjolnir::ExcludedVolumePotential<real_type>::parameter_type>
            >;
        using lj_interaction  = mjolnir::GlobalPairInteraction<
                traits_type, mjolnir::LennardJonesPotential<real_type>,
                mjolnir::NaivePairCalculation<traits_type,
                    typename mjolnir::LennardJonesPotential<real_type>::parameter_type>
            >;

        std::map<std::type_index, bool> found;
        found[typeid(exv_interaction)] = false;
        found[typeid(lj_interaction)]  = false;

        for(const auto& interaction_ptr : ff)
        {
            BOOST_TEST(static_cast<bool>(interaction_ptr));
            const auto& deref = *interaction_ptr;
            found[typeid(deref)] = true;
        }
        BOOST_TEST(found[typeid(exv_interaction)]);
        BOOST_TEST(found[typeid(lj_interaction)]);
    }
}
