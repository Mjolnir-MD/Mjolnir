#define BOOST_TEST_MODULE "test_read_global_stoichiometric_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_global_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_global_stoichiometric_uniform_cubic_pan)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_read_global_stoichiometric_uniform_cubic_pan.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::StoichiometricUniformCubicPanPotential<real_type>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Stoichiometric"
            potential   = "StoichiometricUniformCubicPan"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            epsilon = 3.14
            v0      = 4.0
            range   = 2.0
            particle_kinds = [
            {name = "A", coef = 3},
            {name = "B", coef = 2}
            ]
            parameters = [
                {index = 0, name = "A"},
                {index = 1, name = "B"}
            ]
        )"_toml;

        const auto base = mjolnir::read_stoichiometric_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::GlobalStoichiometricInteraction<traits_type, potential_type>*
            >(base.get());
        BOOST_TEST(static_cast<bool>(derv));
    }

}
