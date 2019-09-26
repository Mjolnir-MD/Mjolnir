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

    using bead_kind = mjolnir::parameter_3SPN2::bead_kind;

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using potential_type      = mjolnir::ThreeSPN2ExcludedVolumePotential<real_type>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "Pair"
            potential = "3SPN2ExcludedVolume"
            ignore.particles_within.bond            = 1
            ignore.particles_within.3SPN2_angle     = 1
            ignore.particles_within.3SPN2_dihedral  = 1
            ignore.particles_within.next_nucleotide = 1
            spatial_partition = {type = "VerletList", margin = 1.0}
            parameters = [
                {index = 0, kind = "P"},
                {index = 1, kind = "S"},
                {index = 2, kind = "A"},
                {index = 3, kind = "T"},
                {index = 4, kind = "G"},
                {index = 5, kind = "C"},
            ]
        )"_toml;

        const auto base = mjolnir::read_global_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<
            mjolnir::GlobalPairInteraction<traits_type, potential_type>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));

        BOOST_TEST_REQUIRE(derv->potential().participants().size() == 6u);

        BOOST_TEST(derv->potential().parameters().at(0) == bead_kind::Phosphate);
        BOOST_TEST(derv->potential().parameters().at(1) == bead_kind::Sugar    );
        BOOST_TEST(derv->potential().parameters().at(2) == bead_kind::BaseA    );
        BOOST_TEST(derv->potential().parameters().at(3) == bead_kind::BaseT    );
        BOOST_TEST(derv->potential().parameters().at(4) == bead_kind::BaseG    );
        BOOST_TEST(derv->potential().parameters().at(5) == bead_kind::BaseC    );
    }
}
