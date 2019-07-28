#define BOOST_TEST_MODULE "test_read_3spn2_base_base_interaction"

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
    using real_type   = typename traits_type::real_type;
    using potential_type      = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;
    using pair_parameter_type = typename potential_type::pair_parameter_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "3SPN2BaseBase"
            spatial_partition = {type = "VerletList", margin = 1.0}
            parameters = [
            {nucleotide_index = 0, S = 0, B =  1, base = "A", B5 = "none", B3 =      4},
            {nucleotide_index = 1, S = 3, B =  4, base = "T", B5 =      1, B3 =      7},
            {nucleotide_index = 2, S = 6, B =  7, base = "G", B5 =      4, B3 =     10},
            {nucleotide_index = 3, S = 9, B = 10, base = "C", B5 =      7, B3 = "none"},
            ]
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::ThreeSPN2BaseBaseInteraction<
                traits_type, mjolnir::VerletList<traits_type, pair_parameter_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        BOOST_TEST_REQUIRE(derv->potential().participants().size() == 4u);

        BOOST_TEST(derv->potential().parameters().at( 1).strand_index == 0u);
        BOOST_TEST(derv->potential().parameters().at( 4).strand_index == 1u);
        BOOST_TEST(derv->potential().parameters().at( 7).strand_index == 2u);
        BOOST_TEST(derv->potential().parameters().at(10).strand_index == 3u);

        BOOST_TEST(derv->potential().parameters().at( 1).S_idx == 0u);
        BOOST_TEST(derv->potential().parameters().at( 4).S_idx == 3u);
        BOOST_TEST(derv->potential().parameters().at( 7).S_idx == 6u);
        BOOST_TEST(derv->potential().parameters().at(10).S_idx == 9u);

        BOOST_TEST(derv->potential().parameters().at( 1).B3_idx ==  4u);
        BOOST_TEST(derv->potential().parameters().at( 4).B3_idx ==  7u);
        BOOST_TEST(derv->potential().parameters().at( 7).B3_idx == 10u);
        BOOST_TEST(derv->potential().parameters().at(10).B3_idx == potential_type::invalid());

        BOOST_TEST(derv->potential().parameters().at( 1).B5_idx == potential_type::invalid());
        BOOST_TEST(derv->potential().parameters().at( 4).B5_idx == 1u);
        BOOST_TEST(derv->potential().parameters().at( 7).B5_idx == 4u);
        BOOST_TEST(derv->potential().parameters().at(10).B5_idx == 7u);

        BOOST_TEST(derv->potential().parameters().at( 1).base == base_kind::A);
        BOOST_TEST(derv->potential().parameters().at( 4).base == base_kind::T);
        BOOST_TEST(derv->potential().parameters().at( 7).base == base_kind::G);
        BOOST_TEST(derv->potential().parameters().at(10).base == base_kind::C);
    }
}
