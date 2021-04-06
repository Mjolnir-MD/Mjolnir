#define BOOST_TEST_MODULE "test_global_stoichiometric_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteraction.hpp>
#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteractionPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(GlobalStoichiometricInteraction_one_on_one_double)
{
    mjolnir::LoggerManager::set_default_logger(
        "test_global_stoichiometric_interaction_one_on_one.log");

    using traits_type   = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type     = traits_type::real_type;
    using coord_type    = traits_type::coordinate_type;
    using boundary_type = traits_type::boundary_type;
    using system_type   = mjolnir::System<traits_type>;
    using topology_type = mjolnir::Topology;

    using potential_type = mjolnir::GlobalStoichiometricInteractionPotential<traits_type>;
    using partition_type = mjolnir::NaivePairCalculation<traits_type, potential_type>;

    using interaction_type = mjolnir::GlobalStoichiometricInteraction<traits_type>;

    // set parameters for test
    constexpr real_type tol = 1e-4;
    constexpr real_type dr  = 1e-5;

    // set parameters for system
    real_type particle_radius   = 1.0;
    real_type interaction_range = 1.0;
    real_type epsilon           = 1.0;

    // generate systems, interactions and potentials.
    std::size_t first_coef  = 1;
    std::size_t second_coef = 1;

    system_type system = system_type(4, boundary_type{});

    system.mass(0)  = 1.0;
    system.mass(1)  = 1.0;
    system.mass(2)  = 1.0;
    system.mass(3)  = 1.0;
    system.rmass(0) = 1.0;
    system.rmass(1) = 1.0;
    system.rmass(2) = 1.0;
    system.rmass(3) = 1.0;
    system.name(0)  = "A";
    system.name(1)  = "B";
    system.name(2)  = "A";
    system.name(3)  = "B";
    system.group(0) = "NONE";
    system.group(1) = "NONE";
    system.group(2) = "NONE";
    system.group(3) = "NONE";

    system.position(0) = coord_type( 0.0, 0.0, 0.0);
    system.position(1) = coord_type( 1.5, 0.0, 0.0);
    system.position(2) = coord_type(-1.5, 0.0, 0.0);
    system.position(3) = coord_type( 3.0, 0.0, 0.0);
    system.velocity(0) = coord_type( 0.0, 0.0, 0.0);
    system.velocity(1) = coord_type( 0.0, 0.0, 0.0);
    system.velocity(2) = coord_type( 0.0, 0.0, 0.0);
    system.velocity(3) = coord_type( 0.0, 0.0, 0.0);
    system.force(0)    = coord_type( 0.0, 0.0, 0.0);
    system.force(1)    = coord_type( 0.0, 0.0, 0.0);
    system.force(2)    = coord_type( 0.0, 0.0, 0.0);
    system.force(3)    = coord_type( 0.0, 0.0, 0.0);

    potential_type potential = potential_type(particle_radius, interaction_range,
        std::vector<std::size_t>{0, 2}, std::vector<std::size_t>{1, 3},
        {}, typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type   ({})
        );

    interaction_type interaction = interaction_type(potential_type{potential},
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()),
        epsilon, first_coef, second_coef);

    topology_type topology = topology_type(4);

    topology.construct_molecules();
    interaction.initialize(system, topology);

    const system_type init = system;
    // check the initial force of the system
    interaction.calc_force(system);
    // check the force between different kind particle
    BOOST_TEST(mjolnir::math::X(system.force(0)) ==  1.5, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::X(system.force(1)) == -1.5, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(system.force(0)) ==  0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(system.force(1)) ==  0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(system.force(0)) ==  0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(system.force(1)) ==  0.0, boost::test_tools::tolerance(tol));

    // check the force between same kind particle is 0
    BOOST_TEST(mjolnir::math::X(system.force(2)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::X(system.force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(system.force(2)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(system.force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(system.force(2)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(system.force(3)) == 0.0, boost::test_tools::tolerance(tol));

    for(std::size_t idx=0; idx<system.size(); ++idx)
    {
        {
            for(std::size_t step=0; step<10; ++step)
            {
                // ----------------------------------------------------------------
                // reset positions
                system = init;
                // move particle for test
                mjolnir::math::X(system.position(idx)) += 0.1 * step;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(system);

                mjolnir::math::X(system.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(system);

                mjolnir::math::X(system.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(system);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::X(system.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
        {
            for(std::size_t step=0; step<10; ++step) {
                // ----------------------------------------------------------------
                // reset positions
                system = init;
                // move particle for test
                mjolnir::math::Y(system.position(idx)) += 0.1 * step;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(system);

                mjolnir::math::Y(system.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(system);
                mjolnir::math::Y(system.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(system);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::Y(system.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
        {
            for(std::size_t step=0; step<10; ++step)
            {
                // ----------------------------------------------------------------
                // reset positions
                system = init;
                // move particle for test
                mjolnir::math::Z(system.position(idx)) += 0.1 * step;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(system);

                mjolnir::math::Z(system.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(system);

                mjolnir::math::Z(system.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(system);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::Z(system.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(GlobalStoichiometricInteraction_total_force_double)
{
    mjolnir::LoggerManager::set_default_logger(
        "test_global_stoichiometric_interaction_total_force_double.log");

    using traits_type   = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type     = traits_type::real_type;
    using coord_type    = traits_type::coordinate_type;
    using boundary_type = traits_type::boundary_type;
    using system_type   = mjolnir::System<traits_type>;
    using topology_type = mjolnir::Topology;

    using potential_type = mjolnir::GlobalStoichiometricInteractionPotential<traits_type>;
    using partition_type = mjolnir::NaivePairCalculation<traits_type, potential_type>;

    using interaction_type = mjolnir::GlobalStoichiometricInteraction<traits_type>;

    // set parameters for test
    constexpr real_type tol = 1e-4;

    // set parameters for system
    real_type particle_radius   = 1.0;
    real_type interaction_range = 1.0;
    real_type epsilon           = 1.0;
    std::size_t particle_a_num  = 5;
    std::size_t particle_b_num  = 5;

    // generate systems, interactions and potentials.
    std::size_t first_coef  = 1;
    std::size_t second_coef = 1;

    system_type system = system_type(particle_a_num+particle_b_num, boundary_type{});

    // generate system configuration
    std::mt19937 mt(123456789);
    real_type box_edge = 2.0;
    std::uniform_real_distribution<real_type> uni(-box_edge, box_edge);

    // particle A setup
    for(std::size_t i = 0; i < particle_a_num; ++i)
    {
        system.mass(i)  = 1.0;
        system.rmass(i) = 1.0;
        system.name(i)  = "A";
        system.group(i) = "NONE";
        system.position(i) = coord_type(uni(mt), uni(mt), uni(mt));
        system.velocity(i) = coord_type(0, 0, 0);
        system.force(i)    = coord_type(0, 0, 0);
    }

    // particle B setup
    for(std::size_t i = particle_a_num; i < particle_a_num+particle_b_num; ++i)
    {
        system.mass(i)  = 1.0;
        system.rmass(i) = 1.0;
        system.name(i)  = "B";
        system.group(i) = "NONE";
        system.position(i) = coord_type(uni(mt), uni(mt), uni(mt));
        system.velocity(i) = coord_type(0, 0, 0);
        system.force(i)    = coord_type(0, 0, 0);

    }

    std::vector<std::size_t> particle_a_arr(particle_a_num);
    std::iota(particle_a_arr.begin(), particle_a_arr.end(), 0);
    std::vector<std::size_t> particle_b_arr(particle_b_num);
    std::iota(particle_b_arr.begin(), particle_b_arr.end(), particle_a_num);
    potential_type potential = potential_type(particle_radius, interaction_range,
        std::move(particle_a_arr), std::move(particle_b_arr),
        {}, typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type   ({})
        );

    interaction_type interaction = interaction_type(potential_type{potential},
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()),
        epsilon, first_coef, second_coef);

    topology_type topology = topology_type(500);

    topology.construct_molecules();
    interaction.initialize(system, topology);

    // check total force is 0
    interaction.calc_force(system);
    coord_type total_force(0.0, 0.0, 0.0);
    for(std::size_t idx=0; idx<system.size(); ++idx)
    {
        total_force += system.force(idx);
    }
    BOOST_TEST(mjolnir::math::X(total_force) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(total_force) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(total_force) == 0.0, boost::test_tools::tolerance(tol));
}

BOOST_AUTO_TEST_CASE(GlobalStoichiometricInteraction_one_on_two_double)
{
    mjolnir::LoggerManager::set_default_logger(
        "test_global_stoichiometric_interaction_one_on_two.log");

    using traits_type   = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type     = traits_type::real_type;
    using coord_type    = traits_type::coordinate_type;
    using boundary_type = traits_type::boundary_type;
    using system_type   = mjolnir::System<traits_type>;
    using topology_type = mjolnir::Topology;

    using potential_type = mjolnir::GlobalStoichiometricInteractionPotential<traits_type>;
    using partition_type = mjolnir::NaivePairCalculation<traits_type, potential_type>;

    using interaction_type = mjolnir::GlobalStoichiometricInteraction<traits_type>;

    // set parameters for test
    constexpr real_type tol = 1e-4;
    constexpr real_type dr  = 1e-5;

    // set parameters for system
    real_type particle_radius   = 1.0;
    real_type interaction_range = 1.0;
    real_type epsilon           = 1.0;

    // generate system, interaction and potential.
    system_type      system = system_type(5, boundary_type{});

    //system for A:B = 1:2 case
    std::size_t first_coef  = 1;
    std::size_t second_coef = 2;

    system.mass(0)  = 1.0;
    system.mass(1)  = 1.0;
    system.mass(2)  = 1.0;
    system.mass(3)  = 1.0;
    system.mass(4)  = 1.0;
    system.rmass(0) = 1.0;
    system.rmass(1) = 1.0;
    system.rmass(2) = 1.0;
    system.rmass(3) = 1.0;
    system.rmass(4) = 1.0;
    system.name(0)  = "A";
    system.name(1)  = "B";
    system.name(2)  = "B";
    system.name(3)  = "A";
    system.name(4)  = "B";
    system.group(0) = "NONE";
    system.group(1) = "NONE";
    system.group(2) = "NONE";
    system.group(3) = "NONE";
    system.group(4) = "NONE";

    system.position(0) = coord_type( 0.0,  0.0, 0.0);
    system.position(1) = coord_type( 1.5,  0.0, 0.0);
    system.position(2) = coord_type( 0.0,  1.5, 0.0);
    system.position(3) = coord_type(-1.0, -1.0, 0.0);
    system.position(4) = coord_type( 3.0,  0.0, 0.0);
    system.velocity(0) = coord_type( 0.0,  0.0, 0.0);
    system.velocity(1) = coord_type( 0.0,  0.0, 0.0);
    system.velocity(2) = coord_type( 0.0,  0.0, 0.0);
    system.velocity(3) = coord_type( 0.0,  0.0, 0.0);
    system.velocity(4) = coord_type( 0.0,  0.0, 0.0);
    system.force(0)    = coord_type( 0.0,  0.0, 0.0);
    system.force(1)    = coord_type( 0.0,  0.0, 0.0);
    system.force(2)    = coord_type( 0.0,  0.0, 0.0);
    system.force(3)    = coord_type( 0.0,  0.0, 0.0);
    system.force(4)    = coord_type( 0.0,  0.0, 0.0);

    potential_type potential = potential_type(particle_radius, interaction_range,
        std::vector<std::size_t>{0, 3}, std::vector<std::size_t>{1, 2, 4},
        {}, typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type   ({})
        );

    interaction_type interaction = interaction_type(potential_type{potential},
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()),
        epsilon, first_coef, second_coef);

    topology_type topology = topology_type(5);

    topology.construct_molecules();
    interaction.initialize(system, topology);

    const system_type init = system;
    // check the force between same kind particle is 0
    interaction.calc_force(system);
    BOOST_TEST(mjolnir::math::X(system.force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::X(system.force(4)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(system.force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(system.force(4)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(system.force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(system.force(4)) == 0.0, boost::test_tools::tolerance(tol));

    for(std::size_t idx=0; idx<system.size(); ++idx)
    {
        {
            for(std::size_t step=0; step<10; ++step)
            {
                // ----------------------------------------------------------------
                // reset positions
                system = init;
                // move particle for test
                mjolnir::math::X(system.position(idx)) += 0.1 * step;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(system);

                mjolnir::math::X(system.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(system);

                mjolnir::math::X(system.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(system);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::X(system.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
        {
            for(std::size_t step=0; step<10; ++step) {
                // ----------------------------------------------------------------
                // reset positions
                system = init;
                // move particle for test
                mjolnir::math::Y(system.position(idx)) += 0.1 * step;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(system);

                mjolnir::math::Y(system.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(system);
                mjolnir::math::Y(system.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(system);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::Y(system.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
        {
            for(std::size_t step=0; step<10; ++step)
            {
                // ----------------------------------------------------------------
                // reset positions
                system = init;
                // move particle for test
                mjolnir::math::Z(system.position(idx)) += 0.1 * step;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(system);

                mjolnir::math::Z(system.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(system);

                mjolnir::math::Z(system.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(system);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE / dr == mjolnir::math::Z(system.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
    }
}
