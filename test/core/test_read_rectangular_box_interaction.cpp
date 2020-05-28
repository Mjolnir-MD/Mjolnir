#define BOOST_TEST_MODULE "test_read_rectangular_box_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_external_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_rectangular_box_interaction_EXV)
{
    mjolnir::LoggerManager::set_default_logger("test_read_rectangular_box_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "RectangularBox"
            potential   = "ExcludedVolumeWall"
            box.lower   = [ 0.0,  0.0,  0.0]
            box.upper   = [10.0, 10.0, 10.0]
            box.margin  = 0.4
            epsilon     = 1.0
            parameters  = [
                {index =   0, radius =  1.0},
                {index = 100, radius = 10.0},
            ]
            )"_toml;

        const auto base = mjolnir::read_external_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::RectangularBoxInteraction<
            traits_type, mjolnir::ExcludedVolumeWallPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        const auto& interaction = *derv;

        BOOST_TEST(interaction.potential().epsilon() == 1.0);

        BOOST_TEST(interaction.potential().participants().size() ==   2u);
        BOOST_TEST(interaction.potential().participants().at(0)  ==   0u);
        BOOST_TEST(interaction.potential().participants().at(1)  == 100u);

        BOOST_TEST(interaction.potential().parameters().at(  0) == real_type( 1.0), boost::test_tools::tolerance(1e-8));
        BOOST_TEST(interaction.potential().parameters().at(100) == real_type(10.0), boost::test_tools::tolerance(1e-8));

        BOOST_TEST(mjolnir::math::X(interaction.lower()) == real_type(0.0), boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Y(interaction.lower()) == real_type(0.0), boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Z(interaction.lower()) == real_type(0.0), boost::test_tools::tolerance(1e-8));

        BOOST_TEST(mjolnir::math::X(interaction.upper()) == real_type(10.0), boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Y(interaction.upper()) == real_type(10.0), boost::test_tools::tolerance(1e-8));
        BOOST_TEST(mjolnir::math::Z(interaction.upper()) == real_type(10.0), boost::test_tools::tolerance(1e-8));
    }
}
