#define BOOST_TEST_MODULE "test_omp_global_pair_debye_huckel_interaction"

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
#include <mjolnir/omp/GlobalPairInteraction.hpp>
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_GlobalPair_DebyeHuckel_calc_force)
{
    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_global_pair_debye_huckel_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::DebyeHuckelPotential<traits_type>;
    using parameter_type   = typename potential_type::parameter_type;
    using partition_type   = mjolnir::UnlimitedGridCellList<traits_type, potential_type>;
    using interaction_type = mjolnir::GlobalPairInteraction<traits_type, potential_type>;
    using rng_type         = mjolnir::RandomNumberGenerator<traits_type>;

    using sequencial_traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using sequencial_potential_type   = mjolnir::DebyeHuckelPotential<sequencial_traits_type>;
    using sequencial_system_type      = mjolnir::System<sequencial_traits_type>;
    using sequencial_partition_type   = mjolnir::UnlimitedGridCellList<sequencial_traits_type, sequencial_potential_type>;
    using sequencial_interaction_type = mjolnir::GlobalPairInteraction<sequencial_traits_type, sequencial_potential_type>;

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

        potential_type potential(potential_type::default_cutoff(), parameters, {},
            typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type   ({}));

        sequencial_potential_type seq_potential(
            potential_type::default_cutoff(), parameters, {},
            typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type   ({}));

        rng_type    rng(123456789);
        system_type sys(N_particle, boundary_type{});
        topology_type topol(N_particle);
        topol.construct_molecules();

        sys.attribute("temperature")    = 300.0;
        sys.attribute("ionic_strength") =   0.2;
        potential.update(sys, topol);

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

        // add perturbation
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            mjolnir::math::X(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Y(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Z(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
        }

        // init sequential one with the same coordinates
        sequencial_system_type seq_sys(N_particle, boundary_type{});
        seq_sys.attribute("temperature")    = 300.0;
        seq_sys.attribute("ionic_strength") =   0.2;
        seq_potential.update(seq_sys, topol);

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

        interaction_type interaction(std::move(potential),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));
        sequencial_interaction_type seq_interaction(std::move(seq_potential),
            mjolnir::SpatialPartition<sequencial_traits_type, sequencial_potential_type>(
                mjolnir::make_unique<sequencial_partition_type>()));

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

    }
}

BOOST_AUTO_TEST_CASE(omp_GlobalPair_DebyeHuckel_calc_force_and_energy)
{
    mjolnir::LoggerManager::set_default_logger("test_omp_global_pair_debye_huckel_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = typename traits_type::real_type;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::DebyeHuckelPotential<traits_type>;
    using parameter_type   = typename potential_type::parameter_type;
    using partition_type   = mjolnir::UnlimitedGridCellList<traits_type, potential_type>;
    using interaction_type = mjolnir::GlobalPairInteraction<traits_type, potential_type>;
    using rng_type         = mjolnir::RandomNumberGenerator<traits_type>;

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

        potential_type potential(potential_type::default_cutoff(), parameters, {},
            typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type   ({}));

        rng_type    rng(123456789);
        system_type sys(N_particle, boundary_type{});
        topology_type topol(N_particle);
        topol.construct_molecules();

        sys.attribute("temperature")    = 300.0;
        sys.attribute("ionic_strength") =   0.2;
        potential.update(sys, topol);

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

        // add perturbation
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            mjolnir::math::X(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Y(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Z(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
        }

        interaction_type interaction(std::move(potential),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));

        interaction.initialize(sys, topol);

        system_type ref_sys = sys;

        constexpr real_type tol = 1e-4;

        // calculate forces with openmp
        sys.preprocess_forces();
        const auto energy = interaction.calc_force_and_energy(sys);
        sys.postprocess_forces();

        ref_sys.preprocess_forces();
        interaction.calc_force(ref_sys);
        ref_sys.postprocess_forces();

        const auto ref_energy = interaction.calc_energy(ref_sys);

        BOOST_TEST(ref_energy == energy, boost::test_tools::tolerance(tol));

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            BOOST_TEST(mjolnir::math::X(sys.force(idx)) == mjolnir::math::X(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(sys.force(idx)) == mjolnir::math::Y(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(sys.force(idx)) == mjolnir::math::Z(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
        }
    }
}
