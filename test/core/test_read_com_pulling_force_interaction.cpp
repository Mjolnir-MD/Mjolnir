#define BOOST_TEST_MODULE "test_read_com_pulling_force_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_external_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_com_pulling_force_interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_read_pulling_force_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "PullingForce"
            parameters  = [
                {indices = [0,1,2],    force = [ 1.0, 1.0, 1.0]},
                {indices = " [0, 3) ", force = [ 2.0, 2.0, 2.0]},
                {indices = ["[0, 3) ", "[4, 6]"], force = [ 3.0, 3.0, 3.0]},
            ]
            )"_toml;

        const auto base = mjolnir::read_external_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::PullingForceInteraction<traits_type>*>(base.get());
        BOOST_TEST(static_cast<bool>(derv));

        const auto& interaction = *derv;

        const auto expected1 = std::vector<std::string>{0, 1, 2};
        const auto expected2 = std::vector<std::string>{0, 1, 2};
        const auto expected3 = std::vector<std::string>{0, 1, 2, 4, 5, 6};

        const auto& actual1 = interaction.parameters().at(0).first;
        const auto& actual2 = interaction.parameters().at(1).first;
        const auto& actual3 = interaction.parameters().at(2).first;

        BOOST_TEST_REQUIRE(actual1.size() == expected1.size());
        BOOST_TEST_REQUIRE(actual2.size() == expected2.size());
        BOOST_TEST_REQUIRE(actual3.size() == expected3.size());

        BOOST_TEST_REQUIRE(std::equal(actual1.begin(), actual1.end(), expected1.begin()));
        BOOST_TEST_REQUIRE(std::equal(actual2.begin(), actual2.end(), expected2.begin()));
        BOOST_TEST_REQUIRE(std::equal(actual3.begin(), actual3.end(), expected3.begin()));

        BOOST_TEST(mjolnir::math::X(interaction.parameters().at(0).second) == 1.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Y(interaction.parameters().at(0).second) == 1.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Z(interaction.parameters().at(0).second) == 1.0, boost::test_tools::tolerance(1e-8));

        BOOST_TEST(mjolnir::math::X(interaction.parameters().at(1).second) == 2.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Y(interaction.parameters().at(1).second) == 2.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Z(interaction.parameters().at(1).second) == 2.0, boost::test_tools::tolerance(1e-8));

        BOOST_TEST(mjolnir::math::X(interaction.parameters().at(2).second) == 3.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Y(interaction.parameters().at(2).second) == 3.0, boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Z(interaction.parameters().at(2).second) == 3.0, boost::test_tools::tolerance(1e-8));
    }
}
