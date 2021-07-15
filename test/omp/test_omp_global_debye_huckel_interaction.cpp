#define BOOST_TEST_MODULE "test_omp_global_pair_debye_huckel_interaction"

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
#include <mjolnir/omp/GlobalPairInteraction.hpp>
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_GlobalPair_DebyeHuckel_calc_force)
{
    namespace test = mjolnir::test;
    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_global_pair_debye_huckel_interaction.log");

    using real_type        = double;
    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::DebyeHuckelPotential<real_type>;
    using parameter_list_type = mjolnir::ParameterList<traits_type, potential_type>;
    using parameter_type   = typename mjolnir::DebyeHuckelParameterList<traits_type>::parameter_type;
    using partition_type   = mjolnir::UnlimitedGridCellList<traits_type, potential_type>;
    using interaction_type = mjolnir::GlobalPairInteraction<traits_type, potential_type>;

    using sequential_traits_type         = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using sequential_parameter_list_type = mjolnir::ParameterList<sequential_traits_type, potential_type>;
    using sequential_parameter_type = typename mjolnir::DebyeHuckelParameterList<sequential_traits_type>::parameter_type;
    using sequential_system_type         = mjolnir::System<sequential_traits_type>;
    using sequential_partition_type      = mjolnir::UnlimitedGridCellList<sequential_traits_type, potential_type>;
    using sequential_interaction_type    = mjolnir::GlobalPairInteraction<sequential_traits_type, potential_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    const std::size_t N_particle = 64;
    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        std::vector<std::pair<std::size_t, parameter_type>> parameters(N_particle);
        for(std::size_t i=0; i<N_particle; ++i)
        {
            parameters[i] = std::make_pair(i, parameter_type{1.0});
        }
        std::vector<std::pair<std::size_t, sequential_parameter_type>> seq_parameters(N_particle);
        for(std::size_t i=0; i<N_particle; ++i)
        {
            seq_parameters[i] = std::make_pair(i, sequential_parameter_type{1.0});
        }

        potential_type potential(5.5);

        mjolnir::DebyeHuckelParameterList<traits_type> rule(parameters, {},
            typename parameter_list_type::ignore_molecule_type("Nothing"),
            typename parameter_list_type::ignore_group_type   ({}));
        parameter_list_type parameter_list(mjolnir::make_unique<
            mjolnir::DebyeHuckelParameterList<traits_type>>(std::move(rule)));

        mjolnir::DebyeHuckelParameterList<sequential_traits_type> seq_rule(seq_parameters, {},
            typename parameter_list_type::ignore_molecule_type("Nothing"),
            typename parameter_list_type::ignore_group_type   ({}));
        sequential_parameter_list_type seq_parameter_list(mjolnir::make_unique<
            mjolnir::DebyeHuckelParameterList<sequential_traits_type>>(std::move(seq_rule)));

        std::mt19937 rng(123456789);
        system_type sys(N_particle, boundary_type{});
        topology_type topol(N_particle);
        topol.construct_molecules();

        sys.attribute("temperature")    = 300.0;
        sys.attribute("ionic_strength") =   0.2;
        potential.initialize(sys);

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto i_x = i % 4;
            const auto i_y = i / 4;
            const auto i_z = i / 16;

            sys.mass(i)     = 1.0;
            sys.position(i) = mjolnir::math::make_coordinate<coordinate_type>(i_x*2.0, i_y*2.0, i_z*2.0);
            sys.velocity(i) = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.force(i)    = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.name(i)     = "X";
            sys.group(i)    = "TEST";
        }

        test::apply_random_perturbation(sys, rng, 0.1);

        // init sequential one with the same coordinates
        sequential_system_type seq_sys(N_particle, boundary_type{});
        seq_sys.attribute("temperature")    = 300.0;
        seq_sys.attribute("ionic_strength") =   0.2;

        assert(sys.size() == seq_sys.size());
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            seq_sys.mass(i)     = sys.mass(i);
            seq_sys.position(i) = sys.position(i);
            seq_sys.velocity(i) = sys.velocity(i);
            seq_sys.force(i)    = sys.force(i);
            seq_sys.name(i)     = sys.name(i);
            seq_sys.group(i)    = sys.group(i);
        }

        interaction_type interaction(potential_type{potential},
            std::move(parameter_list),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));

        sequential_interaction_type seq_interaction(potential_type{potential},
            std::move(seq_parameter_list),
            mjolnir::SpatialPartition<sequential_traits_type, potential_type>(
                mjolnir::make_unique<sequential_partition_type>()));

        interaction    .initialize(sys,     topol);
        seq_interaction.initialize(seq_sys, topol);

        // calculate forces with openmp
        interaction.calc_force(sys);
        sys.postprocess_forces();

        // calculate forces without openmp
        seq_interaction.calc_force(seq_sys);

        // check the values are the same
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            BOOST_TEST(mjolnir::math::X(seq_sys.force(i)) == mjolnir::math::X(sys.force(i)),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(seq_sys.force(i)) == mjolnir::math::Y(sys.force(i)),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(seq_sys.force(i)) == mjolnir::math::Z(sys.force(i)),
                       boost::test_tools::tolerance(tol));
        }
        BOOST_TEST(interaction.calc_energy(sys) == seq_interaction.calc_energy(seq_sys),
                   boost::test_tools::tolerance(tol));

        // check the virials are the same
        for(std::size_t i=0; i<9; ++i)
        {
            BOOST_TEST(sys.virial()[i] == seq_sys.virial()[i], boost::test_tools::tolerance(tol));
        }

        // check calc_force_and_energy
        test::check_force_and_energy(sys, interaction, tol);
    }
}
