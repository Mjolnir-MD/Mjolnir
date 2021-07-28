#define BOOST_TEST_MODULE "test_external_distance_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/external/ExternalDistanceInteraction.hpp>
#include <mjolnir/forcefield/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(ExternalDistanceInteraction_numerical_difference)
{
    namespace test = mjolnir::test;
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
    test::clear_everything(sys);

    std::mt19937 mt(123456789);

    for(int i = 0; i < 10000; ++i)
    {
        test::clear_everything(sys);

        sys.position(0) = coordinate_type(0.0, 0.0, 1.0);
        sys.position(1) = coordinate_type(0.0, 0.0, 1.0);

        test::apply_random_perturbation(sys, mt, 0.01);

        // compare between numerical diff and force implementation
        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, interaction, tol, dr, false);
        test::check_force_and_energy(sys, interaction, tol);
    }
}
