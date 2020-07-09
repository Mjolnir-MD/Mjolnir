#define BOOST_TEST_MODULE "test_external_distance_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/external/ExternalDistanceInteraction.hpp>
#include <mjolnir/forcefield/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(ExternalDistanceInteraction_numerical_difference)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_external_distance_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coordinate_type  = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::ExcludedVolumeWallPotential<real_type>;
    using shape_type       = mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveZDirection<traits_type>>;
    using interaction_type = mjolnir::ExternalDistanceInteraction<traits_type, potential_type, shape_type>;

    interaction_type interaction(shape_type(0.0, 0.5),
        potential_type(/*epsilon = */1.0, /* cutoff = */2.0, {
            {0, 1.0}, {1, 1.0}
        }));

    system_type sys(2, boundary_type{});

    sys.at(0).mass  = 1.0;
    sys.at(1).mass  = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;

    sys.at(0).position = coordinate_type( 0.0, 0.0, 1.0);
    sys.at(1).position = coordinate_type( 0.0, 0.0, 1.0);
    sys.at(0).velocity = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(1).velocity = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(0).force    = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(1).force    = coordinate_type( 0.0, 0.0, 0.0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    std::normal_distribution<real_type> gauss(0.0, 1.0);

    for(int i = 0; i < 10000; ++i)
    {
        sys.at(0).position = coordinate_type(0.0, 0.0, 1.0);
        sys.at(1).position = coordinate_type(0.0, 0.0, 1.0);

        // move particles a bit, randomly. and reset forces.
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.position(idx) += coordinate_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.force(idx)     = coordinate_type(0.0, 0.0, 0.0);
        }
        const system_type init = sys;

        // compare between numerical diff and force implementation
        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::X(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::X(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::X(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::Y(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Y(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::Y(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::Z(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Z(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::Z(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }

    }
}

BOOST_AUTO_TEST_CASE(ExternalDistanceInteraction_calc_force_and_energy)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_external_distance_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coordinate_type  = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::ExcludedVolumeWallPotential<real_type>;
    using shape_type       = mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveZDirection<traits_type>>;
    using interaction_type = mjolnir::ExternalDistanceInteraction<traits_type, potential_type, shape_type>;

    interaction_type interaction(shape_type(0.0, 0.5),
        potential_type(/*epsilon = */1.0, /* cutoff = */2.0, {
            {0, 1.0}, {1, 1.0}
        }));


    system_type sys(2, boundary_type{});

    sys.at(0).mass  = 1.0;
    sys.at(1).mass  = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;

    sys.at(0).position = coordinate_type( 0.0, 0.0, 1.0);
    sys.at(1).position = coordinate_type( 0.0, 0.0, 1.0);
    sys.at(0).velocity = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(1).velocity = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(0).force    = coordinate_type( 0.0, 0.0, 0.0);
    sys.at(1).force    = coordinate_type( 0.0, 0.0, 0.0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    std::normal_distribution<real_type> gauss(0.0, 1.0);

    for(int i = 0; i < 10000; ++i)
    {
        sys.at(0).position = coordinate_type(0.0, 0.0, 1.0);
        sys.at(1).position = coordinate_type(0.0, 0.0, 1.0);

        // move particles a bit, randomly. and reset forces.
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.position(idx) += coordinate_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.force(idx)     = coordinate_type(0.0, 0.0, 0.0);
        }

        system_type ref_sys = sys;

        constexpr real_type tol = 1e-4;

        const auto energy = interaction.calc_force_and_energy(sys);
        const auto ref_energy = interaction.calc_energy(ref_sys);
        interaction.calc_force(ref_sys);
        BOOST_TEST(ref_energy == energy, boost::test_tools::tolerance(tol));

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            BOOST_TEST(mjolnir::math::X(sys.force(idx)) == mjolnir::math::X(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(sys.force(idx)) == mjolnir::math::Y(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(sys.force(idx)) == mjolnir::math::Z(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
        }
    }
}
