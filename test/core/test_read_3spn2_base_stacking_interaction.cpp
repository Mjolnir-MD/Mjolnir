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
            topology    = "none"
            parameters = [
                {S_idx =  0, B5_idx =  1, B3_idx =  4, Base5 = "A", Base3 = "A"},
                {S_idx =  3, B5_idx =  4, B3_idx =  7, Base5 = "A", Base3 = "T"},
                {S_idx =  6, B5_idx =  7, B3_idx = 10, Base5 = "A", Base3 = "G"},
                {S_idx =  9, B5_idx = 10, B3_idx = 13, Base5 = "A", Base3 = "C"},
                {S_idx = 12, B5_idx = 13, B3_idx = 16, Base5 = "T", Base3 = "A"},
                {S_idx = 15, B5_idx = 16, B3_idx = 19, Base5 = "T", Base3 = "T"},
                {S_idx = 18, B5_idx = 19, B3_idx = 22, Base5 = "T", Base3 = "G"},
                {S_idx = 21, B5_idx = 22, B3_idx = 25, Base5 = "T", Base3 = "C"},
                {S_idx = 24, B5_idx = 25, B3_idx = 28, Base5 = "G", Base3 = "A"},
                {S_idx = 27, B5_idx = 28, B3_idx = 31, Base5 = "G", Base3 = "T"},
                {S_idx = 30, B5_idx = 31, B3_idx = 34, Base5 = "G", Base3 = "G"},
                {S_idx = 33, B5_idx = 34, B3_idx = 37, Base5 = "G", Base3 = "C"},
                {S_idx = 36, B5_idx = 37, B3_idx = 40, Base5 = "C", Base3 = "A"},
                {S_idx = 39, B5_idx = 40, B3_idx = 43, Base5 = "C", Base3 = "T"},
                {S_idx = 42, B5_idx = 43, B3_idx = 46, Base5 = "C", Base3 = "G"},
                {S_idx = 45, B5_idx = 46, B3_idx = 49, Base5 = "C", Base3 = "C"},
            ]
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::ThreeSPN2BaseStackingInteraction<traits_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        BOOST_TEST(derv->connection_kind() == "none");
        BOOST_TEST_REQUIRE(derv->parameters().size() == 16u);
        for(std::size_t i=0; i<16; ++i)
        {
            BOOST_TEST(derv->parameters().at(i).first.at(0) == 0 + i*3);
            BOOST_TEST(derv->parameters().at(i).first.at(1) == 1 + i*3);
            BOOST_TEST(derv->parameters().at(i).first.at(2) == 4 + i*3);
        }
        BOOST_TEST(derv->parameters().at( 0).second == base_stack_kind::AA);
        BOOST_TEST(derv->parameters().at( 1).second == base_stack_kind::AT);
        BOOST_TEST(derv->parameters().at( 2).second == base_stack_kind::AG);
        BOOST_TEST(derv->parameters().at( 3).second == base_stack_kind::AC);
        BOOST_TEST(derv->parameters().at( 4).second == base_stack_kind::TA);
        BOOST_TEST(derv->parameters().at( 5).second == base_stack_kind::TT);
        BOOST_TEST(derv->parameters().at( 6).second == base_stack_kind::TG);
        BOOST_TEST(derv->parameters().at( 7).second == base_stack_kind::TC);
        BOOST_TEST(derv->parameters().at( 8).second == base_stack_kind::GA);
        BOOST_TEST(derv->parameters().at( 9).second == base_stack_kind::GT);
        BOOST_TEST(derv->parameters().at(10).second == base_stack_kind::GG);
        BOOST_TEST(derv->parameters().at(11).second == base_stack_kind::GC);
        BOOST_TEST(derv->parameters().at(12).second == base_stack_kind::CA);
        BOOST_TEST(derv->parameters().at(13).second == base_stack_kind::CT);
        BOOST_TEST(derv->parameters().at(14).second == base_stack_kind::CG);
        BOOST_TEST(derv->parameters().at(15).second == base_stack_kind::CC);
    }
}
