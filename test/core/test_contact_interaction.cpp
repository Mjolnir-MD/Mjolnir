#define BOOST_TEST_MODULE "test_contact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/ContactInteraction.hpp>
#include <mjolnir/forcefield/local/GaussianPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(Contact_numerical_difference)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_contact_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::GaussianPotential<real_type>;
    using interaction_type = mjolnir::ContactInteraction<traits_type, potential_type>;

    const real_type k(1.0);
    const real_type native(std::sqrt(3.0));
    const real_type sigma(0.15);

    potential_type   potential(k, sigma, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    std::mt19937 mt(123456789);

    for(std::size_t i = 0; i < 1000; ++i)
    {
        system_type sys(2, boundary_type{});
        test::clear_everything(sys);

        sys.position(0) = coord_type( 0.0, 0.0, 0.0);
        sys.position(1) = coord_type( 1.0, 1.0, 1.0);

        test::apply_random_rotation(sys, mt);
        test::apply_random_perturbation(sys, mt, 0.01);

        constexpr real_type tol = 1e-6;
        constexpr real_type dr  = 1e-6;

        test::check_force(sys, interaction, tol, dr);
        test::check_virial(sys, interaction, tol);
        test::check_force_and_virial(sys, interaction, tol);
        test::check_force_and_energy(sys, interaction, tol);
        test::check_force_energy_virial(sys, interaction, tol);
    }
}
