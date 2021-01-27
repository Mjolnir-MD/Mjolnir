#define BOOST_TEST_MODULE "test_omp_global_stoichiometric_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/math/math.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/UnlimitedGridCellList.hpp>
#include <mjolnir/omp/GlobalStoichiometricInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_GlobalStoichiometric_calc_force)
{
    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_global_stoichiometric_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;

    using potential_type   = mjolnir::GlobalStoichiometricInteractionPotential<traits_type>;
    using partition_type   = mjolnir::UnlimitedGridCellList<traits_type, potential_type>;

    using interaction_type = mjolnir::GlobalStoichiometricInteraction<traits_type>;
    using rng_type         = mjolnir::RandomNumberGenerator<traits_type>;

    using sequencial_traits_type    = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using sequencial_system_type    = mjolnir::System<sequencial_traits_type>;
    using sequencial_potential_type = mjolnir::GlobalStoichiometricInteractionPotential<sequencial_traits_type>;
    using sequencial_partition_type = mjolnir::UnlimitedGridCellList<sequencial_traits_type, sequencial_potential_type>;
    using sequencial_interaction_type = mjolnir::GlobalStoichiometricInteraction<sequencial_traits_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    // paramters for GlobalStoichiometricInteractionPotential
    const std::size_t Na_particle = 32;
    const std::size_t Nb_particle = 32;
    const std::size_t N_partice   = Na_particle + Nb_particle;
    const double      v0          = 2.0;
    const double      range       = 2.0;

    for(int num_thread = 1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        std::vector<std::size_t> a_particle_indices(Na_particle);
        std::iota(a_particle_indices.begin(), a_particle_indices.end(), 0);
        std::vector<std::size_t> b_particle_indices(Nb_particle);
        std::iota(b_particle_indices.begin(), b_particle_indices.end(), Na_particle-1);

        std::vector<std::size_t> seq_a_particle_indices = a_particle_indices;
        std::vector<std::size_t> seq_b_particle_indices = b_particle_indices;

        potential_type potential(
            v0, range, std::move(a_particle_indices), std::move(b_particle_indices), {},
            typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type({}));

        sequencial_potential_type seq_potential(
            v0, range, std::move(seq_a_particle_indices), std::move(seq_b_particle_indices), {},
            typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type({}));

        rng_type    rng(123456789);
        system_type sys(N_partice, boundary_type{});
        topology_type topol(N_partice);

        // locate 64 particle to 9x9x9 size cube
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto i_x = i % 4;
            const auto i_y = i / 4;
            const auto i_z = i / 16;

            sys.mass(i)  = 1.0;
            sys.rmass(i) = 1.0;
            sys.position(i) = mjolnir::math::make_coordinate<coordinate_type>(i_x*3.0, i_y*3.0, i_z*3.0);
            sys.velocity(i) = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.force(i)    = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.name(i)     = "X";
            sys.group(i)    = "TEST";
        }

        // add perturbation
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            mjolnir::math::X(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Y(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Z(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
        }
        potential.update(sys, topol);

        // init sequential one with the same coordinates
        sequencial_system_type seq_sys(N_partice, boundary_type{});
        assert(sys.size() == seq_sys.size());
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            seq_sys.mass(i)     = sys.mass(i);
            seq_sys.rmass(i)    = sys.rmass(i);
            seq_sys.position(i) = sys.position(i);
            seq_sys.velocity(i) = sys.velocity(i);
            seq_sys.force(i)    = sys.force(i);
            seq_sys.name(i)     = sys.name(i);
            seq_sys.group(i)    = sys.group(i);
        }
        seq_potential.update(seq_sys, topol);

        partition_type            celllist;
        sequencial_partition_type seq_celllist;

        topol.construct_molecules();

        // paramters for GlobalStoichiometricInteraction
        const double epsilon = 1.0;
        const std::size_t stoichiometric_coef_a = 1;
        const std::size_t stoichiometric_coef_b = 1;

        interaction_type interaction(std::move(potential),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()),
            epsilon, stoichiometric_coef_a, stoichiometric_coef_b);
        sequencial_interaction_type seq_interaction(std::move(seq_potential),
            mjolnir::SpatialPartition<sequencial_traits_type, sequencial_potential_type>(
                mjolnir::make_unique<sequencial_partition_type>()),
            epsilon, stoichiometric_coef_a, stoichiometric_coef_b);

        interaction    .initialize(sys, topol);
        seq_interaction.initialize(seq_sys, topol);

        // calculate force with openmp
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
    }
}

BOOST_AUTO_TEST_CASE(omp_GlobalStoichiometric_calc_force_and_energy)
{
    mjolnir::LoggerManager::set_default_logger("test_omp_global_stoichiometric_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;

    using potential_type   = mjolnir::GlobalStoichiometricInteractionPotential<traits_type>;
    using partition_type   = mjolnir::UnlimitedGridCellList<traits_type, potential_type>;

    using interaction_type = mjolnir::GlobalStoichiometricInteraction<traits_type>;
    using rng_type         = mjolnir::RandomNumberGenerator<traits_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    const std::size_t Na_particle = 32;
    const std::size_t Nb_particle = 32;
    const std::size_t N_partice   = Na_particle + Nb_particle;
    const double      v0          = 2.0;
    const double      range       = 2.0;

    for(int num_thread = 1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        std::vector<std::size_t> a_particle_indices(Na_particle);
        std::iota(a_particle_indices.begin(), a_particle_indices.end(), 0);
        std::vector<std::size_t> b_particle_indices(Nb_particle);
        std::iota(b_particle_indices.begin(), b_particle_indices.end(), Na_particle-1);

        potential_type potential(
            v0, range, std::move(a_particle_indices), std::move(b_particle_indices), {},
            typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type({}));

        rng_type      rng(123456789);
        system_type   sys(N_partice, boundary_type{});
        topology_type topol(N_partice);

        // locate 64 particle to 9x9x9 size cube
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto i_x = i % 4;
            const auto i_y = i / 4;
            const auto i_z = i / 16;

            sys.mass(i)     = 1.0;
            sys.rmass(i)    = 1.0;
            sys.position(i) = mjolnir::math::make_coordinate<coordinate_type>(i_x*3.0, i_y*3.0, i_z*3.0);
            sys.velocity(i) = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.force(i)    = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.name(i)     = "X";
            sys.group(i)    = "TEST";
        }

        // add perturbation
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            mjolnir::math::X(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Y(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Z(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
        }
        potential.update(sys, topol);

        topol.construct_molecules();

        // paramters for GlobalStoichiometricInteraction
        const double epsilon = 1.0;
        const std::size_t stoichiometric_coef_a = 1;
        const std::size_t stoichiometric_coef_b = 1;

        interaction_type interaction(std::move(potential),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()),
            epsilon, stoichiometric_coef_a, stoichiometric_coef_b);

        interaction.initialize(sys, topol);

        constexpr double tol = 1e-4;
        auto ref_sys = sys;

        // calculate force with openmp
        const auto ref_ene = interaction.calc_energy(ref_sys);

        ref_sys.preprocess_forces();
        interaction.calc_force(ref_sys);
        ref_sys.postprocess_forces();

        sys.preprocess_forces();
        const auto ene = interaction.calc_force_and_energy(sys);
        sys.postprocess_forces();

        BOOST_TEST(ene == ref_ene, boost::test_tools::tolerance(tol));
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            BOOST_TEST(mjolnir::math::X(ref_sys.force(i)) == mjolnir::math::X(sys.force(i)),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(ref_sys.force(i)) == mjolnir::math::Y(sys.force(i)),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(ref_sys.force(i)) == mjolnir::math::Z(sys.force(i)),
                       boost::test_tools::tolerance(tol));
        }
    }
}
