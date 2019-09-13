#define BOOST_TEST_MODULE "test_read_3spn2_base_stacking_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_local_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_3spn2_base_stacking_interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_read_3spn2_base_stacking_interaction.log");

    using base_stack_kind = mjolnir::parameter_3SPN2::base_stack_kind;

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "3SPN2BaseStacking"
            potential   = "3SPN2"
            topology    = "nucleotide"
            parameters = [
                {strand = 0, nucleotide =  0,         S =  0, B =  1, Base = "A"}, # AA
                {strand = 0, nucleotide =  1, P =  2, S =  3, B =  4, Base = "A"}, # AT
                {strand = 0, nucleotide =  2, P =  5, S =  6, B =  7, Base = "T"}, # TG
                {strand = 0, nucleotide =  3, P =  8, S =  9, B = 10, Base = "G"}, # GC
                {strand = 0, nucleotide =  4, P = 11, S = 12, B = 13, Base = "C"}, # CA
                {strand = 0, nucleotide =  5, P = 14, S = 15, B = 16, Base = "A"}, # AT
                {strand = 0, nucleotide =  6, P = 17, S = 18, B = 19, Base = "T"}, # TG
                {strand = 0, nucleotide =  7, P = 20, S = 21, B = 22, Base = "G"}, # GC
                {strand = 0, nucleotide =  8, P = 23, S = 24, B = 25, Base = "C"}, # CA
                {strand = 0, nucleotide =  9, P = 26, S = 27, B = 28, Base = "A"}, # AT
                {strand = 0, nucleotide = 10, P = 29, S = 30, B = 31, Base = "T"}, # TG
                {strand = 0, nucleotide = 11, P = 32, S = 33, B = 34, Base = "G"}, # GC
                {strand = 0, nucleotide = 12, P = 35, S = 36, B = 37, Base = "C"}, # CA
                {strand = 0, nucleotide = 13, P = 38, S = 39, B = 40, Base = "A"}, # AT
                {strand = 0, nucleotide = 14, P = 41, S = 42, B = 43, Base = "T"}, # TG
                {strand = 0, nucleotide = 15, P = 44, S = 45, B = 46, Base = "G"}, # GC
                {strand = 0, nucleotide = 16, P = 47, S = 48, B = 49, Base = "C"}, # x
            ]
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::ThreeSPN2BaseStackingInteraction<traits_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        BOOST_TEST(derv->connection_kind() == "nucleotide");
        BOOST_TEST_REQUIRE(derv->parameters().size() == 16u);
        for(std::size_t i=0; i<16; ++i)
        {
            BOOST_TEST(derv->parameters().at(i).first.at(0) == 0 + i*3);
            BOOST_TEST(derv->parameters().at(i).first.at(1) == 1 + i*3);
            BOOST_TEST(derv->parameters().at(i).first.at(2) == 4 + i*3);
        }
        BOOST_TEST(derv->parameters().at( 0).second == base_stack_kind::AA);
        BOOST_TEST(derv->parameters().at( 1).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at( 2).second == base_stack_kind::TG);
        BOOST_TEST(derv->parameters().at( 3).second == base_stack_kind::GC);
        BOOST_TEST(derv->parameters().at( 4).second == base_stack_kind::CA);
        BOOST_TEST(derv->parameters().at( 5).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at( 6).second == base_stack_kind::TG);
        BOOST_TEST(derv->parameters().at( 7).second == base_stack_kind::GC);
        BOOST_TEST(derv->parameters().at( 8).second == base_stack_kind::CA);
        BOOST_TEST(derv->parameters().at( 9).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at(10).second == base_stack_kind::TG);
        BOOST_TEST(derv->parameters().at(11).second == base_stack_kind::GC);
        BOOST_TEST(derv->parameters().at(12).second == base_stack_kind::CA);
        BOOST_TEST(derv->parameters().at(13).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at(14).second == base_stack_kind::TG);
        BOOST_TEST(derv->parameters().at(15).second == base_stack_kind::GC);
    }
}
