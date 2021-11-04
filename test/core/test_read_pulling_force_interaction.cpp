#define BOOST_TEST_MODULE "test_read_pulling_force_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_external_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_pulling_force_interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_read_pulling_force_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "PullingForce"
            parameters  = [
                {index =   0, force = [ 1.0, 2.0, 10.0]},
                {index = 100, force = [-5.0, 0.0,  0.0]},
                {index = 100, force = 0.0144, direction = [1.0, 2.0, 3.0]},
            ]
            )"_toml;

        const auto base = mjolnir::read_external_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::PullingForceInteraction<traits_type>*>(base.get());
        BOOST_TEST(static_cast<bool>(derv));

        const auto& interaction = *derv;

        BOOST_TEST(interaction.parameters().at(0).first ==   0);
        BOOST_TEST(interaction.parameters().at(1).first == 100);

        BOOST_TEST(mjolnir::math::X(interaction.parameters().at(0).second) ==  1.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Y(interaction.parameters().at(0).second) ==  2.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Z(interaction.parameters().at(0).second) == 10.0, boost::test_tools::tolerance(1e-8));

        BOOST_TEST(mjolnir::math::X(interaction.parameters().at(1).second) ==-5.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Y(interaction.parameters().at(1).second) == 0.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Z(interaction.parameters().at(1).second) == 0.0, boost::test_tools::tolerance(1e-8));

        BOOST_TEST(mjolnir::math::X(interaction.parameters().at(2).second) == 1.0 * 0.0144 / std::sqrt(1.0 + 4.0 + 9.0), boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Y(interaction.parameters().at(2).second) == 2.0 * 0.0144 / std::sqrt(1.0 + 4.0 + 9.0), boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Z(interaction.parameters().at(2).second) == 3.0 * 0.0144 / std::sqrt(1.0 + 4.0 + 9.0), boost::test_tools::tolerance(1e-8));
    }
}
