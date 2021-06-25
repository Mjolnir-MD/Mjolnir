#define BOOST_TEST_MODULE "test_hybrid_forcefield"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/hybrid/HybridForceField.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/forcefield/local/BondLengthInteraction.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(MultipleBasin_2Basin_numerical_difference)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_multiple_basin_forcefield.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    const real_type k(10.0);
    const real_type native1(1.5);
    const real_type native2(2.0);

    mjolnir::LocalForceField<traits_type>      loc1, loc2;
    mjolnir::GlobalForceField<traits_type>     glo1, glo2;
    mjolnir::ExternalForceField<traits_type>   ext1, ext2;
    mjolnir::ConstraintForceField<traits_type> cst1, cst2;

    std::vector<std::pair<std::array<std::size_t, 2>, potential_type>> param1, param2;
    param1.emplace_back(std::array<std::size_t, 2>{{0,1}}, potential_type(k, native1));
    param2.emplace_back(std::array<std::size_t, 2>{{0,1}}, potential_type(k, native2));

    loc1.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param1)));
    loc2.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param2)));

    mjolnir::ForceField<traits_type> ff1(std::move(loc1), std::move(glo1), std::move(ext1), std::move(cst1));
    mjolnir::ForceField<traits_type> ff2(std::move(loc2), std::move(glo2), std::move(ext2), std::move(cst2));

    const real_type lambda = 0.5;
    mjolnir::HybridForceField<traits_type> forcefield(lambda,
        mjolnir::make_unique<mjolnir::ForceField<traits_type>>(ff1),
        mjolnir::make_unique<mjolnir::ForceField<traits_type>>(ff2));

    std::mt19937 rng(123456789);
    for(std::size_t i = 0; i < 1000; ++i)
    {
        system_type sys(2, boundary_type{});
        test::clear_everything(sys);

        sys.at(0).position = coord_type( 0.0, 0.0, 0.0);
        sys.at(1).position = coord_type( 1.0, 1.0, 1.0);

        test::apply_random_perturbation(sys, rng, 0.01);
        test::apply_random_rotation(sys, rng);

        forcefield.initialize(sys); // don't forget this

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, forcefield, tol, dr);

        const auto Etest = forcefield.calc_energy(sys);
        const auto E1    = ff1.calc_energy(sys);
        const auto E2    = ff2.calc_energy(sys);

        BOOST_TEST(Etest == lambda * E1 + (1.0 - lambda) * E2, boost::test_tools::tolerance(tol));
    }
}
