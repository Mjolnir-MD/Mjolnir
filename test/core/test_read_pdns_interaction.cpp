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
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "PDNS"
            potential   = "PDNS"
            spatial_partition.type   = "VerletList"
            spatial_partition.margin = 0.4
            sigma = 1.0
            delta = 0.17453
            cutoff = 10.0
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

        const auto& para = derv->parameters();

        BOOST_TEST_REQUIRE(para.cutoff_ratio() == 10.0);

        const auto nil = mjolnir::ProteinDNANonSpecificParameterList<traits_type>::invalid();

        BOOST_TEST_REQUIRE(para.proteins().size() == 3u);
        BOOST_TEST_REQUIRE(para.proteins().at(0) == 10u);
        BOOST_TEST_REQUIRE(para.proteins().at(1) == 11u);
        BOOST_TEST_REQUIRE(para.proteins().at(2) == 12u);
        BOOST_TEST_REQUIRE(para.parameters().at(para.proteins().at(0)).S3 == nil);
        BOOST_TEST_REQUIRE(para.parameters().at(para.proteins().at(1)).S3 == nil);
        BOOST_TEST_REQUIRE(para.parameters().at(para.proteins().at(2)).S3 == nil);

        BOOST_TEST_REQUIRE(para.dnas().size() == 3u);
        BOOST_TEST_REQUIRE(para.dnas().at(0) == 1u);
        BOOST_TEST_REQUIRE(para.dnas().at(1) == 4u);
        BOOST_TEST_REQUIRE(para.dnas().at(2) == 7u);
        BOOST_TEST_REQUIRE(para.parameters().at(para.dnas().at(0)).S3 == 0u);
        BOOST_TEST_REQUIRE(para.parameters().at(para.dnas().at(1)).S3 == 3u);
        BOOST_TEST_REQUIRE(para.parameters().at(para.dnas().at(2)).S3 == 6u);

        BOOST_TEST_REQUIRE(para.contacts().size()   == 3u);
        BOOST_TEST_REQUIRE(para.contacts().at(0).P  == 10u);
        BOOST_TEST_REQUIRE(para.contacts().at(1).P  == 11u);
        BOOST_TEST_REQUIRE(para.contacts().at(2).P  == 12u);
        BOOST_TEST_REQUIRE(para.contacts().at(0).PN == 9u);
        BOOST_TEST_REQUIRE(para.contacts().at(1).PN == 10u);
        BOOST_TEST_REQUIRE(para.contacts().at(2).PN == 11u);
        BOOST_TEST_REQUIRE(para.contacts().at(0).PC == 11u);
        BOOST_TEST_REQUIRE(para.contacts().at(1).PC == 12u);
        BOOST_TEST_REQUIRE(para.contacts().at(2).PC == 13u);

        BOOST_TEST_REQUIRE(para.contacts().at(0).k == -6.0);
        BOOST_TEST_REQUIRE(para.contacts().at(1).k == -6.0);
        BOOST_TEST_REQUIRE(para.contacts().at(2).k == -6.0);

        BOOST_TEST_REQUIRE(para.contacts().at(0).r0 == 5.1);
        BOOST_TEST_REQUIRE(para.contacts().at(1).r0 == 5.2);
        BOOST_TEST_REQUIRE(para.contacts().at(2).r0 == 5.3);

        BOOST_TEST_REQUIRE(para.contacts().at(0).theta0 == 1.5);
        BOOST_TEST_REQUIRE(para.contacts().at(1).theta0 == 1.6);
        BOOST_TEST_REQUIRE(para.contacts().at(2).theta0 == 1.7);

        BOOST_TEST_REQUIRE(para.contacts().at(0).phi0 == 2.0);
        BOOST_TEST_REQUIRE(para.contacts().at(1).phi0 == 2.1);
        BOOST_TEST_REQUIRE(para.contacts().at(2).phi0 == 2.2);
    }
}
