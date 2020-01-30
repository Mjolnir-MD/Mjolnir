#define BOOST_TEST_MODULE "test_read_afm_fitting_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/forcefield/AFMFit/AFMFitInteraction.hpp>
#include <mjolnir/input/read_external_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_afm_fitting_interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_read_afm_fitting_interaction.log");
    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction = "AFMFlexibleFitting"
        k           = 100.0
        gamma       =   1.0
        pixel_x     =  10.0
        pixel_y     =  10.0
        sigma_x     =   2.0
        sigma_y     =   2.0
        length_x    =   5
        length_y    =   5
        z0          = 0.0
        cutoff      = 5.0
        margin      = 0.5
        image       = [
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.5, 1.0, 0.5,
            0.0, 0.0, 1.0, 2.0, 1.0,
            0.0, 0.0, 0.5, 1.0, 0.5,
        ]
        parameters  = [
        {index = 0, radius = 1.0},
        {index = 1, radius = 2.0},
        {index = 4, radius = 3.0},
        {index = 5, radius = 4.0},
        ]
    )"_toml;

    const auto base = mjolnir::read_external_interaction<traits_type>(v);
    BOOST_TEST(static_cast<bool>(base));

    const auto derv = dynamic_cast<
        const mjolnir::AFMFitInteraction<traits_type>*>(base.get());
    BOOST_TEST(static_cast<bool>(derv));

    BOOST_TEST(derv->k() == 100.0);
    BOOST_TEST(derv->gamma   () ==  1.0);
    BOOST_TEST(derv->pixel_x () == 10.0);
    BOOST_TEST(derv->pixel_y () == 10.0);
    BOOST_TEST(derv->sigma_x () ==  2.0);
    BOOST_TEST(derv->sigma_y () ==  2.0);
    BOOST_TEST(derv->length_x() ==  5u );
    BOOST_TEST(derv->length_y() ==  5u );
    BOOST_TEST(derv->z0      () ==  0.0);
    BOOST_TEST(derv->cutoff  () ==  5.0);
    BOOST_TEST(derv->margin  () ==  0.5);

    BOOST_TEST_REQUIRE(derv->participants().size() == 4u);
    BOOST_TEST_REQUIRE(derv->participants().at(0)  == 0u);
    BOOST_TEST_REQUIRE(derv->participants().at(1)  == 1u);
    BOOST_TEST_REQUIRE(derv->participants().at(2)  == 4u);
    BOOST_TEST_REQUIRE(derv->participants().at(3)  == 5u);

    BOOST_TEST(derv->parameters().at(0) == 1.0);
    BOOST_TEST(derv->parameters().at(1) == 2.0);
    BOOST_TEST(derv->parameters().at(4) == 3.0);
    BOOST_TEST(derv->parameters().at(5) == 4.0);

    BOOST_TEST(derv->image().at( 0) == 0.0);
    BOOST_TEST(derv->image().at( 1) == 0.0);
    BOOST_TEST(derv->image().at( 2) == 0.0);
    BOOST_TEST(derv->image().at( 3) == 0.0);
    BOOST_TEST(derv->image().at( 4) == 0.0);
    BOOST_TEST(derv->image().at( 5) == 0.0);
    BOOST_TEST(derv->image().at( 6) == 0.0);
    BOOST_TEST(derv->image().at( 7) == 0.0);
    BOOST_TEST(derv->image().at( 8) == 0.0);
    BOOST_TEST(derv->image().at( 9) == 0.0);
    BOOST_TEST(derv->image().at(10) == 0.0);
    BOOST_TEST(derv->image().at(11) == 0.0);
    BOOST_TEST(derv->image().at(12) == 0.5);
    BOOST_TEST(derv->image().at(13) == 1.0);
    BOOST_TEST(derv->image().at(14) == 0.5);
    BOOST_TEST(derv->image().at(15) == 0.0);
    BOOST_TEST(derv->image().at(16) == 0.0);
    BOOST_TEST(derv->image().at(17) == 1.0);
    BOOST_TEST(derv->image().at(18) == 2.0);
    BOOST_TEST(derv->image().at(19) == 1.0);
    BOOST_TEST(derv->image().at(20) == 0.0);
    BOOST_TEST(derv->image().at(21) == 0.0);
    BOOST_TEST(derv->image().at(22) == 0.5);
    BOOST_TEST(derv->image().at(23) == 1.0);
    BOOST_TEST(derv->image().at(24) == 0.5);
}
