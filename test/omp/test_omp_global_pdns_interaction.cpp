#define BOOST_TEST_MODULE "test_omp_global_pdns_interaction"

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
#include <mjolnir/omp/ProteinDNANonSpecificInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_PDNS_calc_force)
{
    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_pdns_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::ProteinDNANonSpecificPotential<traits_type>;
    using partition_type   = mjolnir::UnlimitedGridCellList<traits_type, potential_type>;
    using interaction_type = mjolnir::ProteinDNANonSpecificInteraction<traits_type>;
    using rng_type         = mjolnir::RandomNumberGenerator<traits_type>;

    using sequencial_traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using sequencial_potential_type   = mjolnir::ProteinDNANonSpecificPotential<sequencial_traits_type>;
    using sequencial_system_type      = mjolnir::System<sequencial_traits_type>;
    using sequencial_partition_type   = mjolnir::UnlimitedGridCellList<sequencial_traits_type, sequencial_potential_type>;
    using sequencial_interaction_type = mjolnir::ProteinDNANonSpecificInteraction<sequencial_traits_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    // o 0     5 o--o 6  |
    //  \       /        |
    //   o 1   o 7       |
    //  /       \        |
    // o 2     8 o--o 9  |
    //  \       /        |
    //   o 3   o 10      |
    //  /       \        |
    // o 4    11 o--o 12 |

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        potential_type potential(1.0, 0.1, 5.0,
            std::vector<typename potential_type::contact_parameter_type>{
                typename potential_type::contact_parameter_type{1, 0, 2, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0},
                typename potential_type::contact_parameter_type{2, 1, 3, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0},
                typename potential_type::contact_parameter_type{3, 2, 4, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0}
            },
            std::vector<typename potential_type::dna_index_type>{
                typename potential_type::dna_index_type{ 7, 8},
                typename potential_type::dna_index_type{10,11}
            },
            {},
            typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type({}));

        sequencial_potential_type seq_potential(1.0, 0.1, 5.0,
            std::vector<typename sequencial_potential_type::contact_parameter_type>{
                typename sequencial_potential_type::contact_parameter_type{1, 0, 2, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0},
                typename sequencial_potential_type::contact_parameter_type{2, 1, 3, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0},
                typename sequencial_potential_type::contact_parameter_type{3, 2, 4, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0}
            },
            std::vector<typename sequencial_potential_type::dna_index_type>{
                typename sequencial_potential_type::dna_index_type{ 7, 8},
                typename sequencial_potential_type::dna_index_type{10,11}
            },
            {},
            typename sequencial_potential_type::ignore_molecule_type("Nothing"),
            typename sequencial_potential_type::ignore_group_type({}));

        rng_type    rng(123456789);
        system_type sys(13, boundary_type{});
        topology_type topol(13);
        sys.position( 0) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  3.3, 0.0);
        sys.position( 1) = mjolnir::math::make_coordinate<coordinate_type>( 1.9,  1.6, 0.0);
        sys.position( 2) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  0.0, 0.0);
        sys.position( 3) = mjolnir::math::make_coordinate<coordinate_type>( 1.9, -1.6, 0.0);
        sys.position( 4) = mjolnir::math::make_coordinate<coordinate_type>( 0.0, -3.3, 0.0);
        sys.position( 5) = mjolnir::math::make_coordinate<coordinate_type>( 8.8,  3.3, 0.0);
        sys.position( 6) = mjolnir::math::make_coordinate<coordinate_type>(12.0,  3.3, 0.0);
        sys.position( 7) = mjolnir::math::make_coordinate<coordinate_type>( 6.9,  1.6, 0.0);
        sys.position( 8) = mjolnir::math::make_coordinate<coordinate_type>( 8.8,  0.0, 0.0);
        sys.position( 9) = mjolnir::math::make_coordinate<coordinate_type>(12.0,  0.0, 0.0);
        sys.position(10) = mjolnir::math::make_coordinate<coordinate_type>( 6.9, -1.6, 0.0);
        sys.position(11) = mjolnir::math::make_coordinate<coordinate_type>( 8.8, -3.3, 0.0);
        sys.position(12) = mjolnir::math::make_coordinate<coordinate_type>(12.0, -3.3, 0.0);

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.mass(i)     = 1.0;
            sys.velocity(i) = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.force(i)    = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.name(i)     = "X";
            sys.group(i)    = "TEST";
        }

        topol.construct_molecules();

        // add perturbation
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            mjolnir::math::X(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Y(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Z(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
        }
        potential.update(sys, topol);

        // init sequential one with the same coordinates
        sequencial_system_type seq_sys(13, boundary_type{});
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
        seq_potential.update(seq_sys, topol);

        partition_type            celllist;
        sequencial_partition_type seq_celllist;

        topol.construct_molecules();

        interaction_type interaction(std::move(potential),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));
        sequencial_interaction_type seq_interaction(std::move(seq_potential),
            mjolnir::SpatialPartition<sequencial_traits_type, sequencial_potential_type>(
                mjolnir::make_unique<sequencial_partition_type>()));

        interaction    .initialize(sys, topol);
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

BOOST_AUTO_TEST_CASE(omp_PDNS_calc_force_and_energy)
{
    mjolnir::LoggerManager::set_default_logger("test_omp_pdns_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::ProteinDNANonSpecificPotential<traits_type>;
    using partition_type   = mjolnir::UnlimitedGridCellList<traits_type, potential_type>;
    using interaction_type = mjolnir::ProteinDNANonSpecificInteraction<traits_type>;
    using rng_type         = mjolnir::RandomNumberGenerator<traits_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << max_number_of_threads);

    // o 0     5 o--o 6  |
    //  \       /        |
    //   o 1   o 7       |
    //  /       \        |
    // o 2     8 o--o 9  |
    //  \       /        |
    //   o 3   o 10      |
    //  /       \        |
    // o 4    11 o--o 12 |

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        potential_type potential(1.0, 0.1, 5.0,
            std::vector<typename potential_type::contact_parameter_type>{
                typename potential_type::contact_parameter_type{1, 0, 2, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0},
                typename potential_type::contact_parameter_type{2, 1, 3, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0},
                typename potential_type::contact_parameter_type{3, 2, 4, 6.0, 5.0, 3.14 * 0.5, 3.14 * 0.667, 0.0, 0.0}
            },
            std::vector<typename potential_type::dna_index_type>{
                typename potential_type::dna_index_type{ 7, 8},
                typename potential_type::dna_index_type{10,11}
            },
            {},
            typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type({}));

        rng_type    rng(123456789);
        system_type sys(13, boundary_type{});
        topology_type topol(13);
        sys.position( 0) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  3.3, 0.0);
        sys.position( 1) = mjolnir::math::make_coordinate<coordinate_type>( 1.9,  1.6, 0.0);
        sys.position( 2) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  0.0, 0.0);
        sys.position( 3) = mjolnir::math::make_coordinate<coordinate_type>( 1.9, -1.6, 0.0);
        sys.position( 4) = mjolnir::math::make_coordinate<coordinate_type>( 0.0, -3.3, 0.0);
        sys.position( 5) = mjolnir::math::make_coordinate<coordinate_type>( 8.8,  3.3, 0.0);
        sys.position( 6) = mjolnir::math::make_coordinate<coordinate_type>(12.0,  3.3, 0.0);
        sys.position( 7) = mjolnir::math::make_coordinate<coordinate_type>( 6.9,  1.6, 0.0);
        sys.position( 8) = mjolnir::math::make_coordinate<coordinate_type>( 8.8,  0.0, 0.0);
        sys.position( 9) = mjolnir::math::make_coordinate<coordinate_type>(12.0,  0.0, 0.0);
        sys.position(10) = mjolnir::math::make_coordinate<coordinate_type>( 6.9, -1.6, 0.0);
        sys.position(11) = mjolnir::math::make_coordinate<coordinate_type>( 8.8, -3.3, 0.0);
        sys.position(12) = mjolnir::math::make_coordinate<coordinate_type>(12.0, -3.3, 0.0);

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.mass(i)     = 1.0;
            sys.velocity(i) = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.force(i)    = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.name(i)     = "X";
            sys.group(i)    = "TEST";
        }

        topol.construct_molecules();

        // add perturbation
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            mjolnir::math::X(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Y(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Z(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
        }
        potential.update(sys, topol);

        partition_type celllist;
        topol.construct_molecules();

        interaction_type interaction(std::move(potential),
            mjolnir::SpatialPartition<traits_type, potential_type>(
                mjolnir::make_unique<partition_type>()));

        interaction.initialize(sys, topol);

        constexpr double tol = 1e-4;
        auto ref_sys = sys;

        // calculate forces with openmp
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
            BOOST_TEST(mjolnir::math::X(ref_sys.force(i)) == mjolnir::math::X(sys.force(i)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(ref_sys.force(i)) == mjolnir::math::Y(sys.force(i)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(ref_sys.force(i)) == mjolnir::math::Z(sys.force(i)), boost::test_tools::tolerance(tol));
        }
    }
}
