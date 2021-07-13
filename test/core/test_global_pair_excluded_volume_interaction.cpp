#define BOOST_TEST_MODULE "test_global_pair_excluded_volume_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/forcefield/global/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(GlobalPairExcludedVolumeInteraction_numeric_limits)
{
    namespace test = mjolnir::test;
    mjolnir::LoggerManager::set_default_logger("test_global_pair_excluded_volume_interaction.log");
    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = traits_type::real_type;
    using coordinate_type  = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::ExcludedVolumePotential<real_type>;
    using parameter_list_type = mjolnir::ExcludedVolumeParameterList<traits_type>;
    using parameter_type   = typename parameter_list_type::parameter_type;
    using partition_type   = mjolnir::NaivePairCalculation<traits_type, potential_type>;
    using interaction_type = mjolnir::GlobalPairInteraction<traits_type, potential_type>;

    potential_type potential(/*cutoff_ratio = */2.0, /*epsilon = */1.0);

    parameter_list_type  parameter_list(
        std::vector<std::pair<std::size_t, parameter_type>>{
            {0, parameter_type{1.0}}, {1, parameter_type{1.0}}
        }, {}, typename parameter_list_type::ignore_molecule_type("Nothing"),
               typename parameter_list_type::ignore_group_type   ({})
        );

    interaction_type interaction(std::move(potential),
        mjolnir::ParameterList<traits_type, potential_type>(
            mjolnir::make_unique<parameter_list_type>(std::move(parameter_list))),
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()));

    system_type sys(2, boundary_type{});
    test::clear_everything(sys);
    sys.position(0) = coordinate_type(0.0, 0.0, 0.0);
    sys.position(1) = coordinate_type(0.5, 0.0, 0.0);

    topology_type topol(2);
    topol.construct_molecules();

    interaction.initialize(sys, topol);
    // check if the pair of particles are within the cutoff range.
    BOOST_REQUIRE(interaction.calc_energy(sys) != 0.0);

    std::mt19937 rng(123456789);

    const real_type r_min = 0.5;
    const real_type r_max = 4.0;
    const real_type dr = 1e-3;

    const int max_count = (r_max - r_min) / dr;

    for(int i = 0; i < max_count-1; ++i)
    {
        const real_type dist = r_min + i * dr;

        sys.position(0) = coordinate_type(0,0,0);
        sys.position(1) = coordinate_type(dist,0,0);
        sys.force(0)    = coordinate_type(0,0,0);
        sys.force(1)    = coordinate_type(0,0,0);

        test::apply_random_perturbation(sys, rng, 0.1);
        test::apply_random_rotation(sys, rng);

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, interaction, tol, dr);
        test::check_virial(sys, interaction, tol);
        test::check_force_and_virial(sys, interaction, tol);
        test::check_force_and_energy(sys, interaction, tol);
    }
}

