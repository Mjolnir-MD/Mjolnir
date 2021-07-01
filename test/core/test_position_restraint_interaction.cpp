#define BOOST_TEST_MODULE "test_position_restraint_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/external/PositionRestraintInteraction.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(PositionRestraint_numerical_differentiation)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_position_restraint_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coordinate_type  = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::PositionRestraintInteraction<traits_type, potential_type>;

    const potential_type pot1{/*k =*/1.0, /*r0 =*/ 0.0};
    const potential_type pot2{/*k =*/1.0, /*r0 =*/10.0};

    interaction_type interaction(
        std::vector<std::tuple<std::size_t, coordinate_type, potential_type>>{
            std::make_tuple(0, coordinate_type{0.0,0.0,0.0}, pot1),
            std::make_tuple(1, coordinate_type{0.0,0.0,0.0}, pot2)
        });

    std::mt19937 rng(123456789);

    for(std::size_t trial = 0; trial < 1000; ++trial)
    {
        system_type sys(2, boundary_type{});
        test::clear_everything(sys);

        sys.at(0).position = coordinate_type( 0.0,  0.0,  0.0);
        sys.at(1).position = coordinate_type(10.0,  0.0,  0.0);

        interaction.initialize(sys);

        test::apply_random_perturbation(sys, rng, 0.1);

        constexpr real_type tol = 1e-3;
        constexpr real_type dr  = 1e-4;
        test::check_force(sys, interaction, tol, dr);
        test::check_force_and_energy(sys, interaction, tol);
    }
}
