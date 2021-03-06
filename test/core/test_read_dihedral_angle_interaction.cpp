#define BOOST_TEST_MODULE "test_read_dihedral_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/input/read_local_interaction.hpp>


BOOST_AUTO_TEST_CASE(read_dihedral_angle_go_contact)
{
    mjolnir::LoggerManager::set_default_logger("test_read_dihedral_angle_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "DihedralAngle"
            potential   = "ClementiDihedral"
            topology    = "none"
            parameters  = []
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::DihedralAngleInteraction<
            traits_type, mjolnir::ClementiDihedralPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_dihedral_angle_gaussian)
{
    mjolnir::LoggerManager::set_default_logger("test_read_dihedral_angle_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "DihedralAngle"
            potential   = "Gaussian"
            topology    = "none"
            parameters  = []
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::DihedralAngleInteraction<
            traits_type, mjolnir::PeriodicGaussianPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_dihedral_angle_flexible_local)
{
    mjolnir::LoggerManager::set_default_logger("test_read_dihedral_angle_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "DihedralAngle"
            potential   = "FlexibleLocalDihedral"
            topology    = "none"
            parameters  = []
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::DihedralAngleInteraction<
            traits_type, mjolnir::FlexibleLocalDihedralPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_dihedral_angle_cosine)
{
    mjolnir::LoggerManager::set_default_logger("test_read_dihedral_angle_interaction.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = traits_type::real_type;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "DihedralAngle"
            potential   = "Cosine"
            topology    = "none"
            parameters  = []
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::DihedralAngleInteraction<
            traits_type, mjolnir::CosinePotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}
