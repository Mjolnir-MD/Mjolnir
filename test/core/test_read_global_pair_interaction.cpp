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
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "ExcludedVolume"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            epsilon    = 3.14
            parameters = []
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::GlobalPairInteraction<traits_type, potential_type>*
            >(base.get());
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_global_pair_dh)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_pair_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::DebyeHuckelPotential<traits_type>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "DebyeHuckel"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = [
                {index = 0, charge =  1.0},
                {index = 1, charge = -1.0},
            ]
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::GlobalPairInteraction<traits_type, potential_type>*
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

    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "LennardJones"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            parameters = [
                {index = 0, sigma = 2.0, epsilon = 1.5},
                {index = 1, sigma = 2.0, "ε"     = 1.5},
                {index = 2, "σ"   = 2.0, epsilon = 1.5},
                {index = 3, "σ"   = 2.0, "ε"     = 1.5},
            ]
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::GlobalPairInteraction<traits_type, potential_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_global_pair_uni_lj)
{
    mjolnir::LoggerManager::set_default_logger("test_read_global_pair_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::UniformLennardJonesPotential<real_type>;

    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "UniformLennardJones"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            sigma   = 2.0
            epsilon = 1.5
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::GlobalPairInteraction<traits_type, potential_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}
