#define BOOST_TEST_MODULE "test_read_position_restraint_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_external_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_position_restraint_harmonic)
{
    mjolnir::LoggerManager::set_default_logger("test_read_position_restraint_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "PositionRestraint"
            potential   = "Harmonic"
            parameters  = [
                {position = [0.0, 0.0, 0.0], index = 0, k = 10.0, v0 = 0.0},
            ]
            )"_toml;

        const auto base = mjolnir::read_external_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::PositionRestraintInteraction<
            traits_type, mjolnir::HarmonicRestraintPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}
