#define BOOST_TEST_MODULE "test_omp_position_restraint_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/math/math.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/PositionRestraintInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_PositionRestraint)
{
    namespace test = mjolnir::test;

    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_position_restraint_interaction.log");

    using seq_traits_type  = mjolnir::SimulatorTraits      <double, mjolnir::UnlimitedBoundary>;
    using omp_traits_type  = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = typename omp_traits_type::real_type;
    using coordinate_type  = typename omp_traits_type::coordinate_type;
    using boundary_type    = typename omp_traits_type::boundary_type;

    using potential_type   = mjolnir::HarmonicPotential<real_type>;

    using omp_system_type      = mjolnir::System<omp_traits_type>;
    using omp_interaction_type = mjolnir::PositionRestraintInteraction<omp_traits_type, potential_type>;

    using seq_system_type      = mjolnir::System<seq_traits_type>;
    using seq_interaction_type = mjolnir::PositionRestraintInteraction<seq_traits_type, potential_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    std::mt19937 rng(123456789);

    const std::size_t N_particle = 64;
    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        potential_type potential(1.0, 0.0);

        omp_system_type omp_sys(N_particle, boundary_type{});
        seq_system_type seq_sys(N_particle, boundary_type{});
        test::clear_everything(omp_sys);
        test::clear_everything(seq_sys);

        using positions_type = std::vector<std::tuple<std::size_t, coordinate_type, potential_type>>;
        positions_type positions;

        for(std::size_t i=0; i<omp_sys.size(); ++i)
        {
            const auto i_x = i % 4;
            const auto i_y = i / 4;
            const auto i_z = i / 16;

            omp_sys.position(i) = mjolnir::math::make_coordinate<coordinate_type>(i_x*2.0, i_y*2.0, i_z*2.0);
            positions.emplace_back(i, omp_sys.position(i), potential);
        }
        test::apply_random_perturbation(omp_sys, rng, 0.1);

        // init sequential one with the same coordinates
        for(std::size_t i=0; i<seq_sys.size(); ++i)
        {
            seq_sys.position(i) = omp_sys.position(i);
        }

        omp_interaction_type omp_interaction{positions_type(positions)};
        seq_interaction_type seq_interaction{positions_type(positions)};

        omp_interaction.initialize(omp_sys);
        seq_interaction.initialize(seq_sys);

        test::check_force_consistency              (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_and_energy_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_energy_consistency             (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
    }
}
