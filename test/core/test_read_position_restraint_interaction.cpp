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
                {position = [0.0, 0.0, 0.0], index =   0, k = 10.0, v0 =  0.0},
                {position = [1.0,-1.0,10.0], index = 100, k = 3.14, v0 = 10.0},
            ]
            )"_toml;

        const auto base = mjolnir::read_external_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::PositionRestraintInteraction<
            traits_type, mjolnir::HarmonicPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        const auto& interaction = *derv;

        BOOST_TEST(interaction.potentials().size() == 2u);

        BOOST_TEST(std::get<0>(interaction.potentials().at(0)) ==   0u);
        BOOST_TEST(std::get<0>(interaction.potentials().at(1)) == 100u);

        BOOST_TEST(std::get<1>(interaction.potentials().at(0)).at(0) == 0.0);
        BOOST_TEST(std::get<1>(interaction.potentials().at(0)).at(1) == 0.0);
        BOOST_TEST(std::get<1>(interaction.potentials().at(0)).at(2) == 0.0);
        BOOST_TEST(std::get<1>(interaction.potentials().at(1)).at(0) == 1.0);
        BOOST_TEST(std::get<1>(interaction.potentials().at(1)).at(1) ==-1.0);
        BOOST_TEST(std::get<1>(interaction.potentials().at(1)).at(2) ==10.0);


        BOOST_TEST(std::get<2>(interaction.potentials().at(0)).k()  == 10.0);
        BOOST_TEST(std::get<2>(interaction.potentials().at(1)).k()  == 3.14);
        BOOST_TEST(std::get<2>(interaction.potentials().at(0)).v0() ==  0.0);
        BOOST_TEST(std::get<2>(interaction.potentials().at(1)).v0() == 10.0);

    }
}
