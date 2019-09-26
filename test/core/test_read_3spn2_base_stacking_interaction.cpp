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
                {strand = 0, nucleotide = 15, P = 44, S = 45, B = 46, Base = "G"}, # X

                {strand = 1, nucleotide = 16,         S = 47, B = 48, Base = "A"}, # AA
                {strand = 1, nucleotide = 17, P = 49, S = 50, B = 51, Base = "A"}, # AT
                {strand = 1, nucleotide = 18, P = 52, S = 53, B = 54, Base = "T"}, # TG
                {strand = 1, nucleotide = 19, P = 55, S = 56, B = 57, Base = "G"}, # GC
                {strand = 1, nucleotide = 20, P = 58, S = 59, B = 60, Base = "C"}, # CA
                {strand = 1, nucleotide = 21, P = 61, S = 62, B = 63, Base = "A"}, # AT
                {strand = 1, nucleotide = 22, P = 64, S = 65, B = 66, Base = "T"}, # TG
                {strand = 1, nucleotide = 23, P = 67, S = 68, B = 69, Base = "G"}, # GC
                {strand = 1, nucleotide = 24, P = 70, S = 71, B = 72, Base = "C"}, # CA
                {strand = 1, nucleotide = 25, P = 73, S = 74, B = 75, Base = "A"}, # AT
                {strand = 1, nucleotide = 26, P = 76, S = 77, B = 78, Base = "T"}, # TG
                {strand = 1, nucleotide = 27, P = 79, S = 80, B = 81, Base = "G"}, # GC
                {strand = 1, nucleotide = 28, P = 82, S = 83, B = 84, Base = "C"}, # CA
                {strand = 1, nucleotide = 29, P = 85, S = 86, B = 87, Base = "A"}, # AT
                {strand = 1, nucleotide = 30, P = 88, S = 89, B = 90, Base = "T"}, # TG
                {strand = 1, nucleotide = 31, P = 91, S = 92, B = 93, Base = "G"}, # X
            ]
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::ThreeSPN2BaseStackingInteraction<traits_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        BOOST_TEST(derv->connection_kind() == "nucleotide");
        BOOST_TEST_REQUIRE(derv->parameters().size() == 30u);
        for(std::size_t i=0; i<15; ++i)
        {
            BOOST_TEST(derv->parameters().at(i).first.at(0) == 0 + i*3);
            BOOST_TEST(derv->parameters().at(i).first.at(1) == 1 + i*3);
            BOOST_TEST(derv->parameters().at(i).first.at(2) == 4 + i*3);
        }
        for(std::size_t i=15; i<30; ++i)
        {
            BOOST_TEST(derv->parameters().at(i).first.at(0) == 47 + (i-15)*3);
            BOOST_TEST(derv->parameters().at(i).first.at(1) == 48 + (i-15)*3);
            BOOST_TEST(derv->parameters().at(i).first.at(2) == 51 + (i-15)*3);
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

        BOOST_TEST(derv->parameters().at(15).second == base_stack_kind::AA);
        BOOST_TEST(derv->parameters().at(16).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at(17).second == base_stack_kind::TG);
        BOOST_TEST(derv->parameters().at(18).second == base_stack_kind::GC);
        BOOST_TEST(derv->parameters().at(19).second == base_stack_kind::CA);
        BOOST_TEST(derv->parameters().at(20).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at(21).second == base_stack_kind::TG);
        BOOST_TEST(derv->parameters().at(22).second == base_stack_kind::GC);
        BOOST_TEST(derv->parameters().at(23).second == base_stack_kind::CA);
        BOOST_TEST(derv->parameters().at(24).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at(25).second == base_stack_kind::TG);
        BOOST_TEST(derv->parameters().at(26).second == base_stack_kind::GC);
        BOOST_TEST(derv->parameters().at(27).second == base_stack_kind::CA);
        BOOST_TEST(derv->parameters().at(28).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at(29).second == base_stack_kind::TG);
    }
}
