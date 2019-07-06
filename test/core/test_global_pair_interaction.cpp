#define BOOST_TEST_MODULE "test_global_pair_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/interaction/global/GlobalPairInteraction.hpp>
#include <mjolnir/potential/global/IgnoreMolecule.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(GlobalPairInteraction_double)
{
    mjolnir::LoggerManager::set_default_logger("test_global_pair_interaction.log");
    using traits = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    constexpr static traits::real_type tol = 1e-8;

    using real_type        = traits::real_type;
    using coordinate_type  = traits::coordinate_type;
    using boundary_type    = traits::boundary_type;
    using system_type      = mjolnir::System<traits>;
    using potential_type   = mjolnir::LennardJonesPotential<real_type>;
    using parameter_type   = typename potential_type::parameter_type;
    using partition_type   = mjolnir::NaivePairCalculation<traits, parameter_type>;
    using interaction_type = mjolnir::GlobalPairInteraction<traits, potential_type, partition_type>;

    auto normalize = [](const coordinate_type& v){return v / mjolnir::math::length(v);};

    potential_type   potential(std::vector<std::pair<std::size_t, parameter_type>>{
            {0, {/* sigma = */ 1.0, /* epsilon = */1.2}},
            {1, {/* sigma = */ 1.0, /* epsilon = */1.2}}
        }, {}, typename potential_type::ignore_molecule_type("Nothing"),
               typename potential_type::ignore_group_type   ("Nothing")
        );

    interaction_type interaction(potential_type{potential}, partition_type{});

    const real_type eq_dist = 1.0 * std::pow(2.0, 1.0 / 6.0);

    system_type sys(2, boundary_type{});

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

    sys.topology().construct_molecules();

    sys.name(0)  = "X";
    sys.name(1)  = "X";
    sys.group(0) = "NONE";
    sys.group(1) = "NONE";

    interaction.initialize(sys);

    std::mt19937 rng(123456789);
    std::normal_distribution<real_type> gauss(0.0, 1.0);

    const real_type rc = potential.max_cutoff_length();

    const real_type r_min = 0.5;
    const real_type r_max = rc;
    const real_type dr = 1e-3;

    const int max_count = (r_max - r_min) / dr;

    real_type dist = r_min;
    for(int i = 0; i < max_count-1; ++i)
    {
        sys.position(0) = coordinate_type(0,0,0);
        sys.position(1) = coordinate_type(0,0,0);
        sys.force(0)    = coordinate_type(0,0,0);
        sys.force(1)    = coordinate_type(0,0,0);

        sys.position(1) += dist *
            normalize(coordinate_type(gauss(rng), gauss(rng), gauss(rng)));

        const real_type deriv = potential.derivative(0, 1, dist);
        const real_type coef  = std::abs(deriv);

        interaction.calc_force(sys);

        // check force strength
        const real_type force_strength1 = mjolnir::math::length(sys.force(0));
        const real_type force_strength2 = mjolnir::math::length(sys.force(1));
        BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
        BOOST_TEST(coef == force_strength2, boost::test_tools::tolerance(tol));

        // check force direction
        const real_type dir1 = mjolnir::math::dot_product(
            normalize(sys.force(0)), normalize(sys.position(0) - sys.position(1)));
        const real_type dir2 = mjolnir::math::dot_product( 
            normalize(sys.force(1)), normalize(sys.position(1) - sys.position(0)));

        if(dist < eq_dist) // repulsive
        {
            BOOST_TEST(dir1 == 1.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dir2 == 1.0, boost::test_tools::tolerance(tol));
        }
        else // attractive
        {
            BOOST_TEST(dir1 == -1.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dir2 == -1.0, boost::test_tools::tolerance(tol));
        }

        BOOST_TEST(mjolnir::math::length(sys.force(0) + sys.force(1)) == 0.0,
                   boost::test_tools::tolerance(tol));

        dist += dr;
    }
}
