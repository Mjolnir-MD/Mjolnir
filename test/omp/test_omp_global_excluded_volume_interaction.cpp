#define BOOST_TEST_MODULE "test_omp_global_pair_excluded_volume_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/math/math.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/UnlimitedGridCellList.hpp>
#include <mjolnir/omp/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_GlobalPair_ExcludedVolume_calc_force)
{
    namespace test = mjolnir::test;
    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_global_pair_excluded_volume_interaction.log");

    using omp_traits_type = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using seq_traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = double;
    using coordinate_type  = typename omp_traits_type::coordinate_type;
    using boundary_type    = typename omp_traits_type::boundary_type;
    using topology_type    = mjolnir::Topology;

    using potential_type   = mjolnir::ExcludedVolumePotential<real_type>;

    using omp_system_type         = mjolnir::System<omp_traits_type>;
    using omp_parameter_list_type = mjolnir::ParameterList<omp_traits_type, potential_type>;
    using omp_parameter_type      = mjolnir::ExcludedVolumeParameterList<omp_traits_type>::parameter_type;
    using omp_partition_type      = mjolnir::UnlimitedGridCellList<omp_traits_type, potential_type>;
    using omp_interaction_type    = mjolnir::GlobalPairInteraction<omp_traits_type, potential_type>;

    using seq_parameter_list_type = mjolnir::ParameterList<seq_traits_type, potential_type>;
    using seq_parameter_type      = mjolnir::ExcludedVolumeParameterList<seq_traits_type>::parameter_type;
    using seq_system_type         = mjolnir::System<seq_traits_type>;
    using seq_partition_type      = mjolnir::UnlimitedGridCellList<seq_traits_type, potential_type>;
    using seq_interaction_type    = mjolnir::GlobalPairInteraction<seq_traits_type, potential_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    std::mt19937 rng(123456789);

    const std::size_t N_particle = 64;
    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        std::vector<std::pair<std::size_t, omp_parameter_type>> omp_parameters(N_particle);
        for(std::size_t i=0; i<N_particle; ++i)
        {
            omp_parameters[i] = std::make_pair(i, omp_parameter_type{1.0});
        }
        std::vector<std::pair<std::size_t, seq_parameter_type>> seq_parameters(N_particle);
        for(std::size_t i=0; i<N_particle; ++i)
        {
            seq_parameters[i] = std::make_pair(i, seq_parameter_type{1.0});
        }

        potential_type potential(/*cutoff = */2.0, /* epsilon = */1.0);

        mjolnir::ExcludedVolumeParameterList<omp_traits_type> omp_rule(omp_parameters, {},
            typename omp_parameter_list_type::ignore_molecule_type("Nothing"),
            typename omp_parameter_list_type::ignore_group_type   ({})
            );
        omp_parameter_list_type omp_parameter_list(mjolnir::make_unique<
            mjolnir::ExcludedVolumeParameterList<omp_traits_type>>(std::move(omp_rule)));

        mjolnir::ExcludedVolumeParameterList<seq_traits_type> seq_rule(seq_parameters, {},
            typename seq_parameter_list_type::ignore_molecule_type("Nothing"),
            typename seq_parameter_list_type::ignore_group_type   ({})
            );
        seq_parameter_list_type seq_parameter_list(mjolnir::make_unique<
            mjolnir::ExcludedVolumeParameterList<seq_traits_type>>(std::move(seq_rule)));

        topology_type topol(N_particle);
        topol.construct_molecules();

        omp_system_type omp_sys(N_particle, boundary_type{});
        seq_system_type seq_sys(N_particle, boundary_type{});

        for(std::size_t i=0; i<omp_sys.size(); ++i)
        {
            const auto i_x = i % 4;
            const auto i_y = i / 4;
            const auto i_z = i / 16;

            omp_sys.position(i) = mjolnir::math::make_coordinate<coordinate_type>(i_x*2.0, i_y*2.0, i_z*2.0);
        }
        test::apply_random_perturbation(omp_sys, rng, 0.1);

        omp_parameter_list.update(omp_sys, topol, potential);

        // init sequential one with the same coordinates
        for(std::size_t i=0; i<omp_sys.size(); ++i)
        {
            seq_sys.mass(i)     = omp_sys.mass(i);
            seq_sys.position(i) = omp_sys.position(i);
            seq_sys.velocity(i) = omp_sys.velocity(i);
            seq_sys.force(i)    = omp_sys.force(i);
            seq_sys.name(i)     = omp_sys.name(i);
            seq_sys.group(i)    = omp_sys.group(i);
        }
        seq_parameter_list.update(seq_sys, topol, potential);

        omp_interaction_type omp_interaction(potential_type{potential},
            std::move(omp_parameter_list),
            mjolnir::SpatialPartition<omp_traits_type, potential_type>(
                mjolnir::make_unique<omp_partition_type>()));
        seq_interaction_type seq_interaction(potential_type{potential},
            std::move(seq_parameter_list),
            mjolnir::SpatialPartition<seq_traits_type, potential_type>(
                mjolnir::make_unique<seq_partition_type>()));

        omp_interaction.initialize(omp_sys, topol);
        seq_interaction.initialize(seq_sys, topol);

        test::check_force_consistency              (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_and_energy_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_and_virial_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_energy_virial_consistency(omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_energy_consistency             (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
    }
}
