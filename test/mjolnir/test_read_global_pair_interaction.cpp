#define BOOST_TEST_MODULE "test_read_global_pair_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_global_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_global_pair_exv)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_pair_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::ExcludedVolumePotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;
    {
        const toml::table v = toml::table{
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
        };
        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::GlobalPairInteraction<
            traits_type, potential_type,
            mjolnir::NaivePairCalculation<traits_type, parameter_type>
            >*>(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_global_pair_dh)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_pair_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::DebyeHuckelPotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;
    {
        const toml::table v = toml::table{
            {"interaction",       toml::value("Pair")},
            {"potential",         toml::value("DebyeHuckel")},
            {"spatial_partition", toml::value(toml::table{
                        {"type", toml::value("Naive")}
            })},
            {"ignore",            toml::value(toml::table{
                {"molecule",         toml::value("Nothing")},
                {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
            })},
            {"parameters",        toml::value(toml::array{
                toml::table{{"index", 0}, {"charge",  1.0}},
                toml::table{{"index", 1}, {"charge", -1.0}}
            })}
        };
        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::GlobalPairInteraction<
            traits_type, potential_type,
            mjolnir::NaivePairCalculation<traits_type, parameter_type>
            >*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_global_pair_lj)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_pair_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::LennardJonesPotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;

    {
        const toml::table v = toml::table{
            {"interaction",       toml::value("Pair")},
            {"potential",         toml::value("LennardJones")},
            {"spatial_partition", toml::value(toml::table{
                        {"type", toml::value("Naive")}
            })},
            {"ignore",            toml::value(toml::table{
                {"molecule",         toml::value("Nothing")},
                {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
            })},
            {"parameters",        toml::value(toml::array{
                toml::table{{"index", 0}, {"sigma", 2.0}, {"epsilon", 1.5}},
                toml::table{{"index", 1}, {u8"σ",   2.0}, {u8"ε",     1.5}},
                toml::table{{"index", 2}, {u8"σ",   2.0}, {"epsilon", 1.5}},
                toml::table{{"index", 3}, {"sigma", 2.0}, {u8"ε",     1.5}}
            })}
        };
        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::GlobalPairInteraction<
            traits_type, potential_type,
            mjolnir::NaivePairCalculation<traits_type, parameter_type>
            >*>(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_global_pair_uni_lj)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_pair_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::UniformLennardJonesPotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;

    {
        const toml::table v = toml::table{
            {"interaction",       toml::value("Pair")},
            {"potential",         toml::value("UniformLennardJones")},
            {"spatial_partition", toml::value(toml::table{
                        {"type", toml::value("Naive")}
            })},
            {"ignore",            toml::value(toml::table{
                {"molecule",         toml::value("Nothing")},
                {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
            })},
            {"sigma",   2.0},
            {"epsilon", 1.5},
        };
        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::GlobalPairInteraction<
            traits_type, potential_type,
            mjolnir::NaivePairCalculation<traits_type, parameter_type>
            >*>(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}
