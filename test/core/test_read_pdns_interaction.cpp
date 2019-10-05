#define BOOST_TEST_MODULE "test_read_pdns_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_global_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_pdns_interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_read_pdns_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using potential_type = mjolnir::ProteinDNANonSpecificPotential<real_type>;
    using bead_kind = typename potential_type::bead_kind;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "PDNS"
            potential   = "PDNS"
            spatial_partition.type   = "VerletList"
            spatial_partition.margin = 1.0
            sigma = 1.0
            delta = 0.17453
            parameters  = [
            {index = 1, kind = "DNA",     S3 = 0},
            {index = 4, kind = "DNA",     S3 = 3},
            {index = 7, kind = "DNA",     S3 = 6},

            {index = 10, kind = "Protein", PN =  9, PC = 11, k = -6.0, r0 = 5.1, theta0 = 1.5, phi0 = 2.0},
            {index = 11, kind = "Protein", PN = 10, PC = 12, k = -6.0, r0 = 5.2, theta0 = 1.6, phi0 = 2.1},
            {index = 12, kind = "Protein", PN = 11, PC = 13, k = -6.0, r0 = 5.3, theta0 = 1.7, phi0 = 2.2},
            ]
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::ProteinDNANonSpecificInteraction<traits_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        const auto& pot = derv->potential();

        BOOST_TEST_REQUIRE(pot.participants().size() == 6u);

        BOOST_TEST_REQUIRE(pot.participants().at(0) ==  1u);
        BOOST_TEST_REQUIRE(pot.participants().at(1) ==  4u);
        BOOST_TEST_REQUIRE(pot.participants().at(2) ==  7u);
        BOOST_TEST_REQUIRE(pot.participants().at(3) == 10u);
        BOOST_TEST_REQUIRE(pot.participants().at(4) == 11u);
        BOOST_TEST_REQUIRE(pot.participants().at(5) == 12u);

        BOOST_TEST(pot.parameters().at( 1).kind == bead_kind::DNA);
        BOOST_TEST(pot.parameters().at( 4).kind == bead_kind::DNA);
        BOOST_TEST(pot.parameters().at( 7).kind == bead_kind::DNA);
        BOOST_TEST(pot.parameters().at(10).kind == bead_kind::Protein);
        BOOST_TEST(pot.parameters().at(11).kind == bead_kind::Protein);
        BOOST_TEST(pot.parameters().at(12).kind == bead_kind::Protein);

        BOOST_TEST(pot.parameters().at( 1).S3 == 0u);
        BOOST_TEST(pot.parameters().at( 4).S3 == 3u);
        BOOST_TEST(pot.parameters().at( 7).S3 == 6u);

        BOOST_TEST(pot.parameters().at(10).PN ==  9u);
        BOOST_TEST(pot.parameters().at(11).PN == 10u);
        BOOST_TEST(pot.parameters().at(12).PN == 11u);

        BOOST_TEST(pot.parameters().at(10).PC == 11u);
        BOOST_TEST(pot.parameters().at(11).PC == 12u);
        BOOST_TEST(pot.parameters().at(12).PC == 13u);

        const auto pp1 = pot.prepare_params(1, 10);
        const auto pp2 = pot.prepare_params(1, 11);
        const auto pp3 = pot.prepare_params(1, 12);

        BOOST_TEST(pp1.S3 == 0);
        BOOST_TEST(pp2.S3 == 0);
        BOOST_TEST(pp3.S3 == 0);

        BOOST_TEST(pp1.k == -6.0);
        BOOST_TEST(pp2.k == -6.0);
        BOOST_TEST(pp3.k == -6.0);

        BOOST_TEST(pp1.r0 == 5.1);
        BOOST_TEST(pp2.r0 == 5.2);
        BOOST_TEST(pp3.r0 == 5.3);

        BOOST_TEST(pp1.theta0 == 1.5);
        BOOST_TEST(pp2.theta0 == 1.6);
        BOOST_TEST(pp3.theta0 == 1.7);

        BOOST_TEST(pp1.phi0 == 2.0);
        BOOST_TEST(pp2.phi0 == 2.1);
        BOOST_TEST(pp3.phi0 == 2.2);

        const auto pp4 = pot.prepare_params(4, 10);
        const auto pp5 = pot.prepare_params(4, 11);
        const auto pp6 = pot.prepare_params(4, 12);

        BOOST_TEST(pp4.S3 == 3);
        BOOST_TEST(pp5.S3 == 3);
        BOOST_TEST(pp6.S3 == 3);

        BOOST_TEST(pp4.k == -6.0);
        BOOST_TEST(pp5.k == -6.0);
        BOOST_TEST(pp6.k == -6.0);

        BOOST_TEST(pp4.r0 == 5.1);
        BOOST_TEST(pp5.r0 == 5.2);
        BOOST_TEST(pp6.r0 == 5.3);

        BOOST_TEST(pp4.theta0 == 1.5);
        BOOST_TEST(pp5.theta0 == 1.6);
        BOOST_TEST(pp6.theta0 == 1.7);

        BOOST_TEST(pp4.phi0 == 2.0);
        BOOST_TEST(pp5.phi0 == 2.1);
        BOOST_TEST(pp6.phi0 == 2.2);

    }
}
