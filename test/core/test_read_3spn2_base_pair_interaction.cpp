#define BOOST_TEST_MODULE "test_read_3spn2_base_pair_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_global_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_3spn2_base_base_interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_read_3spn2_base_base_interaction.log");

    using base_kind = mjolnir::parameter_3SPN2::base_kind;

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "3SPN2BasePair"
            potential   = "3SPN2"
            spatial_partition = {type = "VerletList", margin = 1.0}
            parameters = [
            {strand = 0,          S =   0, B =   1, Base = "A"},
            {strand = 0, P =   2, S =   3, B =   4, Base = "T"},
            {strand = 0, P =   5, S =   6, B =   7, Base = "G"},
            {strand = 0, P =   8, S =   9, B =  10, Base = "C"},
            ]
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::ThreeSPN2BasePairInteraction<traits_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        BOOST_TEST_REQUIRE(derv->parameters().participants().size() == 4u);

        BOOST_TEST_REQUIRE(derv->parameters().participants().at(0) ==  1u);
        BOOST_TEST_REQUIRE(derv->parameters().participants().at(1) ==  4u);
        BOOST_TEST_REQUIRE(derv->parameters().participants().at(2) ==  7u);
        BOOST_TEST_REQUIRE(derv->parameters().participants().at(3) == 10u);

        BOOST_TEST(derv->parameters().parameters().at( 1).S_idx == 0u);
        BOOST_TEST(derv->parameters().parameters().at( 4).S_idx == 3u);
        BOOST_TEST(derv->parameters().parameters().at( 7).S_idx == 6u);
        BOOST_TEST(derv->parameters().parameters().at(10).S_idx == 9u);

        BOOST_TEST(derv->parameters().parameters().at( 1).base == base_kind::A);
        BOOST_TEST(derv->parameters().parameters().at( 4).base == base_kind::T);
        BOOST_TEST(derv->parameters().parameters().at( 7).base == base_kind::G);
        BOOST_TEST(derv->parameters().parameters().at(10).base == base_kind::C);
    }
}
