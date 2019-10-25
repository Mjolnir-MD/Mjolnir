#define BOOST_TEST_MODULE "test_read_pwmcos_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_global_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_pwmcos_interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_read_pdns_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "PWMcos"
            potential   = "PWMcos"
            spatial_partition.type   = "VerletList"
            spatial_partition.margin = 0.4
            sigma        = 1.0
            phi          = 0.17453
            energy_unit  = 0.593
            energy_shift = 0.0
            cutoff       = 5.0
            parameters   = [
            {index = 4, kind = "DNA", S = 3, B5 = 1, B3 =  7, base = "A"},
            {index = 7, kind = "DNA", S = 6, B5 = 4, B3 = 10, base = "T"},

            {index = 5, kind = "Protein", CaN =  4, CaC = 6, gamma = -1.0, epsilon = 0.0, r0 = 5.0, theta1_0 = 2.0944, theta2_0 = 3.1415, theta3_0 = 0.7854, A = 0.1, C = 0.2, G = 0.3, T = 0.4},
            {index = 6, kind = "Protein", CaN =  5, CaC = 7, gamma = -1.1, epsilon = 1.0, r0 = 6.0, theta1_0 = 3.1415, theta2_0 = 2.0944, theta3_0 = 2.3562, A = 0.4, C = 0.3, G = 0.2, T = 0.1},
            ]
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST_REQUIRE(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::PWMcosInteraction<traits_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST_REQUIRE(static_cast<bool>(derv));

        const auto& pot = derv->potential();

        BOOST_TEST(pot.cutoff_ratio() == 5.0);
        BOOST_TEST(pot.energy_unit()  == 0.593);
        BOOST_TEST(pot.energy_shift() == 0.0);
        BOOST_TEST(pot.sigma()        == 1.0);
        BOOST_TEST(pot.phi()          == 0.17453);

        const auto nil = mjolnir::PWMcosPotential<traits_type>::invalid();

        const auto A = mjolnir::PWMcosPotential<traits_type>::base_kind::A;
        const auto T = mjolnir::PWMcosPotential<traits_type>::base_kind::T;
        const auto C = mjolnir::PWMcosPotential<traits_type>::base_kind::C;
        const auto G = mjolnir::PWMcosPotential<traits_type>::base_kind::G;

        BOOST_TEST(pot.proteins().size() == 2u);
        BOOST_TEST(pot.proteins().at(0) == 5u);
        BOOST_TEST(pot.proteins().at(1) == 6u);
        BOOST_TEST(pot.parameters().at(pot.proteins().at(0)).dna_index == nil);
        BOOST_TEST(pot.parameters().at(pot.proteins().at(1)).dna_index == nil);

        BOOST_TEST(pot.dnas().size() == 2u);
        BOOST_TEST(pot.dnas().at(0) == 4u);
        BOOST_TEST(pot.dnas().at(1) == 7u);
        BOOST_TEST(pot.parameters().at(pot.dnas().at(0)).dna_index == 0u);
        BOOST_TEST(pot.parameters().at(pot.dnas().at(1)).dna_index == 1u);

        BOOST_TEST(pot.contacts().size()    == 2u);
        BOOST_TEST(pot.contacts().at(0).Ca  == 5u);
        BOOST_TEST(pot.contacts().at(1).Ca  == 6u);
        BOOST_TEST(pot.contacts().at(0).CaN == 4u);
        BOOST_TEST(pot.contacts().at(1).CaN == 5u);
        BOOST_TEST(pot.contacts().at(0).CaC == 6u);
        BOOST_TEST(pot.contacts().at(1).CaC == 7u);

        BOOST_TEST(pot.contacts().at(0).gamma == -1.0);
        BOOST_TEST(pot.contacts().at(1).gamma == -1.1);
        BOOST_TEST(pot.contacts().at(0).epsilon == 0.0);
        BOOST_TEST(pot.contacts().at(1).epsilon == 1.0);

        BOOST_TEST(pot.contacts().at(0).r0 == 5.0);
        BOOST_TEST(pot.contacts().at(1).r0 == 6.0);

        BOOST_TEST(pot.contacts().at(0).theta1_0 == 2.0944);
        BOOST_TEST(pot.contacts().at(1).theta1_0 == 3.1415);
        BOOST_TEST(pot.contacts().at(0).theta2_0 == 3.1415);
        BOOST_TEST(pot.contacts().at(1).theta2_0 == 2.0944);
        BOOST_TEST(pot.contacts().at(0).theta3_0 == 0.7854);
        BOOST_TEST(pot.contacts().at(1).theta3_0 == 2.3562);

        BOOST_TEST(pot.contacts().at(0).PWM[static_cast<std::size_t>(A)] == 0.1);
        BOOST_TEST(pot.contacts().at(0).PWM[static_cast<std::size_t>(C)] == 0.2);
        BOOST_TEST(pot.contacts().at(0).PWM[static_cast<std::size_t>(G)] == 0.3);
        BOOST_TEST(pot.contacts().at(0).PWM[static_cast<std::size_t>(T)] == 0.4);
        BOOST_TEST(pot.contacts().at(1).PWM[static_cast<std::size_t>(A)] == 0.4);
        BOOST_TEST(pot.contacts().at(1).PWM[static_cast<std::size_t>(C)] == 0.3);
        BOOST_TEST(pot.contacts().at(1).PWM[static_cast<std::size_t>(G)] == 0.2);
        BOOST_TEST(pot.contacts().at(1).PWM[static_cast<std::size_t>(T)] == 0.1);

        const auto dna4_pro5 = pot.prepare_params(5, 4);
        const auto dna4_pro6 = pot.prepare_params(6, 4);
        const auto dna7_pro5 = pot.prepare_params(5, 7);
        const auto dna7_pro6 = pot.prepare_params(6, 7);

        BOOST_TEST(dna4_pro5.S == 3u);
        BOOST_TEST(dna4_pro6.S == 3u);
        BOOST_TEST(dna7_pro5.S == 6u);
        BOOST_TEST(dna7_pro6.S == 6u);

        BOOST_TEST(dna4_pro5.B5 == 1u);
        BOOST_TEST(dna4_pro6.B5 == 1u);
        BOOST_TEST(dna7_pro5.B5 == 4u);
        BOOST_TEST(dna7_pro6.B5 == 4u);

        BOOST_TEST(dna4_pro5.B3 == 7u);
        BOOST_TEST(dna4_pro6.B3 == 7u);
        BOOST_TEST(dna7_pro5.B3 == 10u);
        BOOST_TEST(dna7_pro6.B3 == 10u);

        BOOST_TEST(dna4_pro5.base == A);
        BOOST_TEST(dna4_pro6.base == A);
        BOOST_TEST(dna7_pro5.base == T);
        BOOST_TEST(dna7_pro6.base == T);
    }
}
