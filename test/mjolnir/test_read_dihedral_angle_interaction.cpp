#define BOOST_TEST_MODULE "test_read_dihedral_angle_interaction"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/input/read_local_interaction.hpp>
#include <test/util/traits.hpp>

BOOST_AUTO_TEST_CASE(read_dihedral_angle_harmonic)
{
    mjolnir::LoggerManager::set_default_logger("test_read_dihedral_angle_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::test::traits<real_type>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{
            {"interaction", toml::value("DihedralAngle")},
            {"potential",   toml::value("Harmonic")},
            {"topology",    toml::value("bond")},
            {"parameters",  toml::value(toml::array(/*empty*/))}
        };
        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::DihedralAngleInteraction<
            traits_type, mjolnir::HarmonicPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_dihedral_angle_go_contact)
{
    mjolnir::LoggerManager::set_default_logger("test_read_dihedral_angle_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::test::traits<real_type>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{
            {"interaction", toml::value("DihedralAngle")},
            {"potential",   toml::value("ClementiDihedral")},
            {"topology",    toml::value("contact")},
            {"parameters",  toml::value(toml::array(/*empty*/))}
        };
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

    using real_type = double;
    using traits_type = mjolnir::test::traits<real_type>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{
            {"interaction", toml::value("DihedralAngle")},
            {"potential",   toml::value("Gaussian")},
            {"topology",    toml::value("none")},
            {"parameters",  toml::value(toml::array(/*empty*/))}
        };
        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::DihedralAngleInteraction<
            traits_type, mjolnir::AngularGaussianPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}

BOOST_AUTO_TEST_CASE(read_dihedral_angle_flexible_local)
{
    mjolnir::LoggerManager::set_default_logger("test_read_dihedral_angle_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::test::traits<real_type>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{
            {"interaction", toml::value("DihedralAngle")},
            {"potential",   toml::value("FlexibleLocalDihedral")},
            {"topology",    toml::value("none")},
            {"parameters",  toml::value(toml::array(/*empty*/))}
        };
        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv = dynamic_cast<mjolnir::DihedralAngleInteraction<
            traits_type, mjolnir::FlexibleLocalDihedralPotential<real_type>>*
            >(base.get()); // check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
    }
}
