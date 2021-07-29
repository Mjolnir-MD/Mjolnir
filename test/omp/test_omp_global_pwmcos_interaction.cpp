#define BOOST_TEST_MODULE "test_omp_global_pwmcos_interaction"

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
#include <mjolnir/omp/PWMcosInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_PWMcos_calc_force)
{
    namespace test = mjolnir::test;
    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_pwmcos_interaction.log");

    using omp_traits_type = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using seq_traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = double;
    using coordinate_type  = typename omp_traits_type::coordinate_type;
    using boundary_type    = typename omp_traits_type::boundary_type;
    using topology_type    = mjolnir::Topology;

    using potential_type   = mjolnir::PWMcosPotential<real_type>;

    using omp_system_type         = mjolnir::System<omp_traits_type>;
    using omp_parameter_list_type = mjolnir::PWMcosParameterList<omp_traits_type>;
    using omp_interaction_type    = mjolnir::PWMcosInteraction<omp_traits_type>;
    using omp_partition_type      = mjolnir::UnlimitedGridCellList<omp_traits_type, potential_type>;

    using seq_system_type         = mjolnir::System<seq_traits_type>;
    using seq_parameter_list_type = mjolnir::PWMcosParameterList<seq_traits_type>;
    using seq_partition_type      = mjolnir::UnlimitedGridCellList<seq_traits_type, potential_type>;
    using seq_interaction_type    = mjolnir::PWMcosInteraction<seq_traits_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    //         5  6
    //         o--o   |
    //             \ 7|
    // o 0          o |
    //  \      8  9/  |
    //   o 1   o--o   |
    //  /          \10|
    // o 2          o |
    //  \     11 12/  |
    //   o 3   o--o   |
    //  /          \13|
    // o 4          o |
    //        14 15/  |
    //         o--o   |

    std::mt19937 rng(123456789);

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        omp_parameter_list_type omp_parameter_list(1.0, 0.1, 1.0, 0.0, 5.0,
            std::vector<typename omp_parameter_list_type::contact_parameter_type>{
                typename omp_parameter_list_type::contact_parameter_type{1, 0, 2, 6.0, 0.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 3.14 * 0.667, 0.0, {{1.0, 1.0, 1.0, 1.0}} },
                typename omp_parameter_list_type::contact_parameter_type{2, 1, 3, 6.0, 0.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 3.14 * 0.667, 0.0, {{1.0, 1.0, 1.0, 1.0}} },
                typename omp_parameter_list_type::contact_parameter_type{3, 2, 4, 6.0, 0.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 3.14 * 0.667, 0.0, {{1.0, 1.0, 1.0, 1.0}} }
            },
            std::vector<typename omp_parameter_list_type::dna_parameter_type>{
                typename omp_parameter_list_type::dna_parameter_type{omp_parameter_list_type::base_kind::A,  8,  9,  5, 11},
                typename omp_parameter_list_type::dna_parameter_type{omp_parameter_list_type::base_kind::A, 11, 15,  8, 14},
            },
            {},
            typename omp_parameter_list_type::ignore_molecule_type("Nothing"),
            typename omp_parameter_list_type::ignore_group_type({}));

        seq_parameter_list_type seq_parameter_list(1.0, 0.1, 1.0, 0.0, 5.0,
            std::vector<typename seq_parameter_list_type::contact_parameter_type>{
                typename seq_parameter_list_type::contact_parameter_type{1, 0, 2, 6.0, 0.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 3.14 * 0.667, 0.0, {{1.0, 1.0, 1.0, 1.0}} },
                typename seq_parameter_list_type::contact_parameter_type{2, 1, 3, 6.0, 0.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 3.14 * 0.667, 0.0, {{1.0, 1.0, 1.0, 1.0}} },
                typename seq_parameter_list_type::contact_parameter_type{3, 2, 4, 6.0, 0.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 3.14 * 0.667, 0.0, {{1.0, 1.0, 1.0, 1.0}} }
            },
            std::vector<typename seq_parameter_list_type::dna_parameter_type>{
                typename seq_parameter_list_type::dna_parameter_type{seq_parameter_list_type::base_kind::A,  8,  9,  5, 11},
                typename seq_parameter_list_type::dna_parameter_type{seq_parameter_list_type::base_kind::A, 11, 15,  8, 14},
            },
            {},
            typename seq_parameter_list_type::ignore_molecule_type("Nothing"),
            typename seq_parameter_list_type::ignore_group_type({}));

        omp_system_type omp_sys(16, boundary_type{});
        test::clear_everything(omp_sys);

        omp_sys.position( 0) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  3.3, 0.0);
        omp_sys.position( 1) = mjolnir::math::make_coordinate<coordinate_type>( 1.9,  1.6, 0.0);
        omp_sys.position( 2) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  0.0, 0.0);
        omp_sys.position( 3) = mjolnir::math::make_coordinate<coordinate_type>( 1.9, -1.6, 0.0);
        omp_sys.position( 4) = mjolnir::math::make_coordinate<coordinate_type>( 0.0, -3.3, 0.0);
        omp_sys.position( 5) = mjolnir::math::make_coordinate<coordinate_type>( 6.9,  5.0, 0.0);
        omp_sys.position( 6) = mjolnir::math::make_coordinate<coordinate_type>( 8.8,  5.0, 0.0);
        omp_sys.position( 7) = mjolnir::math::make_coordinate<coordinate_type>(12.0,  3.3, 0.0);
        omp_sys.position( 8) = mjolnir::math::make_coordinate<coordinate_type>( 6.9,  1.6, 0.0);
        omp_sys.position( 9) = mjolnir::math::make_coordinate<coordinate_type>( 8.8,  1.6, 0.0);
        omp_sys.position(10) = mjolnir::math::make_coordinate<coordinate_type>(12.0,  0.0, 0.0);
        omp_sys.position(11) = mjolnir::math::make_coordinate<coordinate_type>( 6.9, -1.6, 0.0);
        omp_sys.position(12) = mjolnir::math::make_coordinate<coordinate_type>( 8.8, -1.6, 0.0);
        omp_sys.position(13) = mjolnir::math::make_coordinate<coordinate_type>(12.0, -3.3, 0.0);
        omp_sys.position(14) = mjolnir::math::make_coordinate<coordinate_type>( 6.9, -5.0, 0.0);
        omp_sys.position(15) = mjolnir::math::make_coordinate<coordinate_type>( 8.8, -5.0, 0.0);

        test::apply_random_perturbation(omp_sys, rng, 0.1);

        topology_type topol(16);
        topol.construct_molecules();

        // init seq one with the same coordinates
        seq_system_type seq_sys(16, boundary_type{});
        test::clear_everything(seq_sys);
        for(std::size_t i=0; i<seq_sys.size(); ++i)
        {
            seq_sys.mass(i)     = omp_sys.mass(i);
            seq_sys.position(i) = omp_sys.position(i);
            seq_sys.velocity(i) = omp_sys.velocity(i);
            seq_sys.force(i)    = omp_sys.force(i);
            seq_sys.name(i)     = omp_sys.name(i);
            seq_sys.group(i)    = omp_sys.group(i);
        }

        omp_interaction_type omp_interaction(potential_type{}, std::move(omp_parameter_list),
            mjolnir::SpatialPartition<omp_traits_type, potential_type>(
                mjolnir::make_unique<omp_partition_type>()));

        seq_interaction_type seq_interaction(potential_type{}, std::move(seq_parameter_list),
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
