#define BOOST_TEST_MODULE "test_recutangular_box_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/external/RectangularBoxInteraction.hpp>
#include <mjolnir/forcefield/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(PositionRestraint_Harmonic)
{
    namespace test = mjolnir::test;
    mjolnir::LoggerManager::set_default_logger(
            "test_recutangular_box_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coordinate_type  = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::ExcludedVolumeWallPotential<real_type>;
    using interaction_type = mjolnir::RectangularBoxInteraction<traits_type, potential_type>;

    coordinate_type lower( 0.0,  0.0,  0.0);
    coordinate_type upper(10.0, 10.0, 10.0);

    interaction_type interaction(lower, upper, /*margin = */0.5,
        potential_type(/*epsilon = */1.0, /* cutoff = */2.0, {
            {0, 1.0}, {1, 1.0}
        }));

    system_type sys(2, boundary_type{});

    std::mt19937 mt(123456789);

    for(int i = 0; i < 10000; ++i)
    {
        test::clear_everything(sys);
        sys.at(0).position = coordinate_type( 1.0,  9.0,  1.0);
        sys.at(1).position = coordinate_type( 9.0,  1.0,  9.0);

        test::apply_random_perturbation(sys, mt, 0.01);

        // compare between numerical diff and force implementation
        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, interaction, tol, dr, false);
        test::check_force_and_energy(sys, interaction, tol);
    }
}
