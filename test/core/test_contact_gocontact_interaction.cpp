#define BOOST_TEST_MODULE "test_bond_length_gocontact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/GoContactInteraction.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

// a test for specialization of ContactInteraction for GoContactPotential.

BOOST_AUTO_TEST_CASE(ContactGoContact_numerical_difference)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_contact_gocontact.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::GoContactPotential<real_type>;
    using interaction_type = mjolnir::ContactInteraction<traits_type, potential_type>;

    const real_type k(1.0);
    const real_type native(std::sqrt(3.0));

    potential_type   potential(k, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    std::mt19937 mt(123456789);

    for(std::size_t i = 0; i < 1000; ++i)
    {
        system_type sys(2, boundary_type{});
        test::clear_everything(sys);

        sys.at(0).position = coord_type( 0.0, 0.0, 0.0);
        sys.at(1).position = coord_type( 1.0, 1.0, 1.0);

        test::apply_random_rotation(sys, mt);
        test::apply_random_perturbation(sys, mt, 0.01);

        constexpr real_type tol = 2e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, interaction, tol, dr);
        test::check_virial(sys, interaction, tol);
        test::check_force_and_energy(sys, interaction, tol);
    }
}

BOOST_AUTO_TEST_CASE(ContactGoContactInteraction)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_contact_gocontact.log");
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::GoContactPotential<real_type>;
    using interaction_type = mjolnir::ContactInteraction<traits_type, potential_type>;

    constexpr real_type tol = 1e-7;

    const coord_type pos1(1.0, 0.0, 0.0);
    const coord_type pos2(0.0, 0.0, 0.0);
    system_type sys(2, boundary_type{});
    test::clear_everything(sys);

    const real_type k(1.0);
    const real_type native(1.0);

    potential_type   potential(k, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    const int N = 10000;
    const real_type dx = 0.001;
    for(int i = 1; i < N; ++i)
    {
        const auto dist = 0.6 + i * dx;
        sys.position(0) = coord_type(0.0,  0.0, 0.0);
        sys.position(1) = coord_type(dist, 0.0, 0.0);
        sys.force(0)    = coord_type(0.0, 0.0, 0.0);
        sys.force(1)    = coord_type(0.0, 0.0, 0.0);

        const real_type deriv = potential.derivative(dist);

        interaction.initialize(sys);
        interaction.calc_force(sys);

        const real_type force_strength1 = mjolnir::math::length(sys.force(0));
        const real_type force_strength2 = mjolnir::math::length(sys.force(1));

        BOOST_TEST(force_strength1 == std::abs(deriv), boost::test_tools::tolerance(tol));
        BOOST_TEST(force_strength2 == std::abs(deriv), boost::test_tools::tolerance(tol));

        const auto E_pot = potential.potential(dist);
        const auto E_int = interaction.calc_energy(sys);

        BOOST_TEST(E_pot == E_int, boost::test_tools::tolerance(tol));
    }
}

