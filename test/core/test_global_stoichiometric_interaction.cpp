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

BOOST_AUTO_TEST_CASE(GlobalStoichiometricInteraction_double)
{
    mjolnir::LoggerManager::set_default_logger(
        "test_global_stoichiometric_interaction.log");

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
    std::vector<system_type>      systems;
    std::vector<interaction_type> interactions;
    std::vector<potential_type>   potentials;
    std::vector<topology_type>    topologies;

    // system for A:B = 1:1 case.
    std::size_t first_coef  = 1;
    std::size_t second_coef = 1;

    systems.emplace_back(system_type(4, boundary_type{}));

    systems[0].mass(0)  = 1.0;
    systems[0].mass(1)  = 1.0;
    systems[0].mass(2)  = 1.0;
    systems[0].mass(3)  = 1.0;
    systems[0].rmass(0) = 1.0;
    systems[0].rmass(1) = 1.0;
    systems[0].rmass(2) = 1.0;
    systems[0].rmass(3) = 1.0;
    systems[0].name(0)  = "A";
    systems[0].name(1)  = "B";
    systems[0].name(2)  = "A";
    systems[0].name(3)  = "B";
    systems[0].group(0) = "NONE";
    systems[0].group(1) = "NONE";
    systems[0].group(2) = "NONE";
    systems[0].group(3) = "NONE";

    systems[0].position(0) = coord_type( 0.0, 0.0, 0.0);
    systems[0].position(1) = coord_type( 1.5, 0.0, 0.0);
    systems[0].position(2) = coord_type(-1.5, 0.0, 0.0);
    systems[0].position(3) = coord_type( 3.0, 0.0, 0.0);
    systems[0].velocity(0) = coord_type( 0.0, 0.0, 0.0);
    systems[0].velocity(1) = coord_type( 0.0, 0.0, 0.0);
    systems[0].velocity(2) = coord_type( 0.0, 0.0, 0.0);
    systems[0].velocity(3) = coord_type( 0.0, 0.0, 0.0);
    systems[0].force(0)    = coord_type( 0.0, 0.0, 0.0);
    systems[0].force(1)    = coord_type( 0.0, 0.0, 0.0);
    systems[0].force(2)    = coord_type( 0.0, 0.0, 0.0);
    systems[0].force(3)    = coord_type( 0.0, 0.0, 0.0);

    potentials.emplace_back(potential_type(particle_radius, interaction_range,
        std::vector<std::size_t>{0, 2}, std::vector<std::size_t>{1, 3},
        {}, typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type   ({})
        ));

    interactions.emplace_back(interaction_type(potential_type{potentials[0]},
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()),
        epsilon, first_coef, second_coef));

    topologies.emplace_back(topology_type(4));

    // system for A:B = 1:2 case
    second_coef = 2;
    systems.emplace_back(system_type(5, boundary_type{}));

    systems[1].mass(0)  = 1.0;
    systems[1].mass(1)  = 1.0;
    systems[1].mass(2)  = 1.0;
    systems[1].mass(3)  = 1.0;
    systems[1].mass(4)  = 1.0;
    systems[1].rmass(0) = 1.0;
    systems[1].rmass(1) = 1.0;
    systems[1].rmass(2) = 1.0;
    systems[1].rmass(3) = 1.0;
    systems[1].rmass(4) = 1.0;
    systems[1].name(0)  = "A";
    systems[1].name(1)  = "B";
    systems[1].name(2)  = "B";
    systems[1].name(3)  = "A";
    systems[1].name(4)  = "B";
    systems[1].group(0) = "NONE";
    systems[1].group(1) = "NONE";
    systems[1].group(2) = "NONE";
    systems[1].group(3) = "NONE";
    systems[1].group(4) = "NONE";

    systems[1].position(0) = coord_type( 0.0,  0.0, 0.0);
    systems[1].position(1) = coord_type( 1.5,  0.0, 0.0);
    systems[1].position(2) = coord_type( 0.0,  1.5, 0.0);
    systems[1].position(3) = coord_type(-1.0, -1.0, 0.0);
    systems[1].position(4) = coord_type( 3.0,  0.0, 0.0);
    systems[1].velocity(0) = coord_type( 0.0,  0.0, 0.0);
    systems[1].velocity(1) = coord_type( 0.0,  0.0, 0.0);
    systems[1].velocity(2) = coord_type( 0.0,  0.0, 0.0);
    systems[1].velocity(3) = coord_type( 0.0,  0.0, 0.0);
    systems[1].velocity(4) = coord_type( 0.0,  0.0, 0.0);
    systems[1].force(0)    = coord_type( 0.0,  0.0, 0.0);
    systems[1].force(1)    = coord_type( 0.0,  0.0, 0.0);
    systems[1].force(2)    = coord_type( 0.0,  0.0, 0.0);
    systems[1].force(3)    = coord_type( 0.0,  0.0, 0.0);
    systems[1].force(4)    = coord_type( 0.0,  0.0, 0.0);

    potentials.emplace_back(potential_type(particle_radius, interaction_range,
        std::vector<std::size_t>{0, 3}, std::vector<std::size_t>{1, 2, 4},
        {}, typename potential_type::ignore_molecule_type("Nothing"),
            typename potential_type::ignore_group_type   ({})
        ));

    interactions.emplace_back(interaction_type(potential_type{potentials[1]},
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()),
        epsilon, first_coef, second_coef));

    topologies.emplace_back(topology_type(5));

    for(std::size_t system_idx=0; system_idx<systems.size(); ++system_idx)
    {
        topologies[system_idx].construct_molecules();
        interactions[system_idx].initialize(systems[system_idx], topologies[system_idx]);
    }

    // check the force between same kind particle is 0
    system_type sys0 = systems[0];
    interactions[0].calc_force(sys0);
    BOOST_TEST(mjolnir::math::X(systems[0].force(2)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::X(systems[0].force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(systems[0].force(2)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(systems[0].force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(systems[0].force(2)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(systems[0].force(3)) == 0.0, boost::test_tools::tolerance(tol));

    system_type sys1 = systems[1];
    interactions[1].calc_force(sys1);
    BOOST_TEST(mjolnir::math::X(systems[1].force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::X(systems[1].force(4)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(systems[1].force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(systems[1].force(4)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(systems[1].force(3)) == 0.0, boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(systems[1].force(4)) == 0.0, boost::test_tools::tolerance(tol));

    for(std::size_t system_idx=0; system_idx<1; ++system_idx)
    {
        system_type&      sys         = systems[system_idx];
        const system_type init        = sys;
        interaction_type& interaction = interactions[system_idx];
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            {
                for(std::size_t step=0; step<10; ++step)
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;
                    // move particle for test
                    mjolnir::math::X(sys.position(idx)) += 0.1 * step;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::X(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::X(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::X(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
            }
            {
                for(std::size_t step=0; step<10; ++step) {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;
                    // move particle for test
                    mjolnir::math::Y(sys.position(idx)) += 0.1 * step;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::Y(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);
                    mjolnir::math::Y(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::Y(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
            }
            {
                for(std::size_t step=0; step<10; ++step)
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;
                    // move particle for test
                    mjolnir::math::Z(sys.position(idx)) += 0.1 * step;

                    // calc U(x-dx)
                    const auto E0 = interaction.calc_energy(sys);

                    mjolnir::math::Z(sys.position(idx)) += dr;

                    // calc F(x)
                    interaction.calc_force(sys);

                    mjolnir::math::Z(sys.position(idx)) += dr;

                    // calc U(x+dx)
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE / dr == mjolnir::math::Z(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
            }
        }
    }
}
