#define BOOST_TEST_MODULE "test_omp_external_distance_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/math/math.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/forcefield/external/LennardJonesWallPotential.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/UnlimitedGridCellList.hpp>
#include <mjolnir/omp/ExternalDistanceInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_ExternalDistacne_calc_force)
{
    namespace test = mjolnir::test;

    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_external_distance_interaction.log");

    using seq_traits_type  = mjolnir::SimulatorTraits      <double, mjolnir::UnlimitedBoundary>;
    using omp_traits_type  = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = typename omp_traits_type::real_type;
    using coordinate_type  = typename omp_traits_type::coordinate_type;
    using boundary_type    = typename omp_traits_type::boundary_type;
    using topology_type    = mjolnir::Topology;

    using potential_type   = mjolnir::LennardJonesWallPotential<real_type>;
    using parameter_type   = typename potential_type::parameter_type;

    using omp_system_type      = mjolnir::System<omp_traits_type>;
    using omp_shape_type       = mjolnir::AxisAlignedPlane<omp_traits_type, mjolnir::PositiveZDirection<omp_traits_type>>;
    using omp_interaction_type = mjolnir::ExternalDistanceInteraction<omp_traits_type, potential_type, omp_shape_type>;

    using seq_system_type      = mjolnir::System<seq_traits_type>;
    using seq_shape_type       = mjolnir::AxisAlignedPlane<seq_traits_type, mjolnir::PositiveZDirection<seq_traits_type>>;
    using seq_interaction_type = mjolnir::ExternalDistanceInteraction<seq_traits_type, potential_type, seq_shape_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    const std::size_t N_particle = 64;
    std::vector<std::pair<std::size_t, parameter_type>> parameters(N_particle);
    for(std::size_t i=0; i<N_particle; ++i)
    {
        parameters.emplace_back(i, parameter_type(1.0, 1.0));
    }
    const potential_type potential(2.5, parameters);
    topology_type topol(N_particle);
    topol.construct_molecules();

    std::mt19937 rng(123456789);

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        omp_system_type omp_sys(N_particle, boundary_type{});
        seq_system_type seq_sys(N_particle, boundary_type{});
        test::clear_everything(omp_sys);
        test::clear_everything(seq_sys);

        for(std::size_t i=0; i<omp_sys.size(); ++i)
        {
            const auto i_x = i % 4;
            const auto i_y = i / 4;
            const auto i_z = i / 16;

            omp_sys.position(i) = mjolnir::math::make_coordinate<coordinate_type>(i_x*2.0, i_y*2.0, i_z*2.0);
        }

        test::apply_random_perturbation(omp_sys, rng, 0.1);

        // init sequential one with the same coordinates
        for(std::size_t i=0; i<seq_sys.size(); ++i)
        {
            seq_sys.mass(i)     = omp_sys.mass(i);
            seq_sys.position(i) = omp_sys.position(i);
            seq_sys.velocity(i) = omp_sys.velocity(i);
            seq_sys.force(i)    = omp_sys.force(i);
            seq_sys.name(i)     = omp_sys.name(i);
            seq_sys.group(i)    = omp_sys.group(i);
        }

        omp_shape_type omp_xyplane(0.0);
        seq_shape_type seq_xyplane(0.0);

        omp_xyplane.initialize(omp_sys, potential);
        seq_xyplane.initialize(seq_sys, potential);

        omp_interaction_type omp_interaction(std::move(omp_xyplane), potential_type(potential));
        seq_interaction_type seq_interaction(std::move(seq_xyplane), potential_type(potential));

        omp_interaction.initialize(omp_sys);
        seq_interaction.initialize(seq_sys);

        test::check_force_consistency              (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_and_energy_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_energy_consistency             (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
    }
}
