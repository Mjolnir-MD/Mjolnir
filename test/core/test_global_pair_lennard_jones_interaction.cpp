#define BOOST_TEST_MODULE "test_global_pair_lennard_jones_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/forcefield/global/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(GlobalPairLennardJonesInteraction_numeric_limits)
{
    mjolnir::LoggerManager::set_default_logger("test_global_pair_lennard_jones_interaction.log");
    using traits = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = traits::real_type;
    using coordinate_type  = traits::coordinate_type;
    using boundary_type    = traits::boundary_type;
    using system_type      = mjolnir::System<traits>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::LennardJonesPotential<traits>;
    using parameter_type   = typename potential_type::parameter_type;
    using partition_type   = mjolnir::NaivePairCalculation<traits, potential_type>;
    using interaction_type = mjolnir::GlobalPairInteraction<traits, potential_type>;

    auto normalize = [](const coordinate_type& v){return v / mjolnir::math::length(v);};

    potential_type   potential(potential_type::default_cutoff(),
        std::vector<std::pair<std::size_t, parameter_type>>{
            {0, {/* sigma = */ 1.0, /* epsilon = */1.2}},
            {1, {/* sigma = */ 1.0, /* epsilon = */1.2}}
        }, {}, typename potential_type::ignore_molecule_type("Nothing"),
               typename potential_type::ignore_group_type   ({})
        );

    interaction_type interaction(potential_type{potential},
        mjolnir::SpatialPartition<traits, potential_type>(
            mjolnir::make_unique<partition_type>()));

    system_type sys(2, boundary_type{});
    topology_type topol(2);

    sys.mass(0)  = 1.0;
    sys.mass(1)  = 1.0;
    sys.rmass(0) = 1.0;
    sys.rmass(1) = 1.0;

    sys.position(0) = coordinate_type(0.0, 0.0, 0.0);
    sys.position(1) = coordinate_type(0.5, 0.0, 0.0);
    sys.velocity(0) = coordinate_type(0.0, 0.0, 0.0);
    sys.velocity(1) = coordinate_type(0.0, 0.0, 0.0);
    sys.force(0)    = coordinate_type(0.0, 0.0, 0.0);
    sys.force(1)    = coordinate_type(0.0, 0.0, 0.0);

    topol.construct_molecules();

    sys.name(0)  = "X";
    sys.name(1)  = "X";
    sys.group(0) = "NONE";
    sys.group(1) = "NONE";

    interaction.initialize(sys, topol);

    std::mt19937 rng(123456789);
    std::normal_distribution<real_type> gauss(0.0, 1.0);

    const real_type rc = potential.max_cutoff_length();

    const real_type r_min = 0.5;
    const real_type r_max = rc;
    const real_type dr = 1e-3;

    const int max_count = (r_max - r_min) / dr;

    for(int i = 0; i < max_count-1; ++i)
    {
        const real_type dist = r_min + i * dr;

        sys.position(0) = coordinate_type(0,0,0);
        sys.position(1) = coordinate_type(0,0,0);
        sys.force(0)    = coordinate_type(0,0,0);
        sys.force(1)    = coordinate_type(0,0,0);

        sys.position(1) += dist *
            normalize(coordinate_type(gauss(rng), gauss(rng), gauss(rng)));

        const auto init = sys;

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

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
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

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
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

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

        // -----------------------------------------------------------------
        // check virial

        sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.force(idx) = coord_type(0,0,0);
        }
        interaction.calc_force(sys);

        matrix33_type vir(0,0,0, 0,0,0, 0,0,0);
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            vir += math::tensor_product(sys.position(idx), sys.force(idx));
        }

        BOOST_TEST(sys.virial()(0,0) == vir(0,0), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(0,1) == vir(0,1), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(0,2) == vir(0,2), boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.virial()(1,0) == vir(1,0), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(1,1) == vir(1,1), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(1,2) == vir(1,2), boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.virial()(2,0) == vir(2,0), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(2,1) == vir(2,1), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(2,2) == vir(2,2), boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(GlobalPairLennardJonesInteraction_numeric_force_and_energy)
{
    mjolnir::LoggerManager::set_default_logger("test_global_pair_lennard_jones_interaction.log");
    using traits = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = traits::real_type;
    using coordinate_type  = traits::coordinate_type;
    using boundary_type    = traits::boundary_type;
    using system_type      = mjolnir::System<traits>;
    using topology_type    = mjolnir::Topology;
    using potential_type   = mjolnir::LennardJonesPotential<traits>;
    using parameter_type   = typename potential_type::parameter_type;
    using partition_type   = mjolnir::NaivePairCalculation<traits, potential_type>;
    using interaction_type = mjolnir::GlobalPairInteraction<traits, potential_type>;

    auto normalize = [](const coordinate_type& v){return v / mjolnir::math::length(v);};

    potential_type   potential(potential_type::default_cutoff(),
        std::vector<std::pair<std::size_t, parameter_type>>{
            {0, {/* sigma = */ 1.0, /* epsilon = */1.2}},
            {1, {/* sigma = */ 1.0, /* epsilon = */1.2}}
        }, {}, typename potential_type::ignore_molecule_type("Nothing"),
               typename potential_type::ignore_group_type   ({})
        );

    interaction_type interaction(potential_type{potential},
        mjolnir::SpatialPartition<traits, potential_type>(
            mjolnir::make_unique<partition_type>()));

    system_type sys(2, boundary_type{});
    topology_type topol(2);

    sys.mass(0)  = 1.0;
    sys.mass(1)  = 1.0;
    sys.rmass(0) = 1.0;
    sys.rmass(1) = 1.0;

    sys.position(0) = coordinate_type(0.0, 0.0, 0.0);
    sys.position(1) = coordinate_type(0.5, 0.0, 0.0);
    sys.velocity(0) = coordinate_type(0.0, 0.0, 0.0);
    sys.velocity(1) = coordinate_type(0.0, 0.0, 0.0);
    sys.force(0)    = coordinate_type(0.0, 0.0, 0.0);
    sys.force(1)    = coordinate_type(0.0, 0.0, 0.0);

    topol.construct_molecules();

    sys.name(0)  = "X";
    sys.name(1)  = "X";
    sys.group(0) = "NONE";
    sys.group(1) = "NONE";

    interaction.initialize(sys, topol);

    std::mt19937 rng(123456789);
    std::normal_distribution<real_type> gauss(0.0, 1.0);

    const real_type rc = potential.max_cutoff_length();

    const real_type r_min = 0.5;
    const real_type r_max = rc;
    const real_type dr = 1e-3;

    const int max_count = (r_max - r_min) / dr;

    for(int i = 0; i < max_count-1; ++i)
    {
        const real_type dist = r_min + i * dr;

        sys.position(0) = coordinate_type(0,0,0);
        sys.position(1) = coordinate_type(0,0,0);
        sys.force(0)    = coordinate_type(0,0,0);
        sys.force(1)    = coordinate_type(0,0,0);

        sys.position(1) += dist *
            normalize(coordinate_type(gauss(rng), gauss(rng), gauss(rng)));

        system_type ref_sys = sys;

        constexpr real_type tol = 1e-4;

        const auto energy = interaction.calc_force_and_energy(sys);
        const auto ref_energy = interaction.calc_energy(ref_sys);
        interaction.calc_force(ref_sys);
        BOOST_TEST(ref_energy == energy, boost::test_tools::tolerance(tol));

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            BOOST_TEST(mjolnir::math::X(sys.force(idx)) == mjolnir::math::X(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(sys.force(idx)) == mjolnir::math::Y(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(sys.force(idx)) == mjolnir::math::Z(ref_sys.force(idx)), boost::test_tools::tolerance(tol));
        }
    }
}
