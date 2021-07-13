#define BOOST_TEST_MODULE "test_PDNS_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificInteraction.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>

BOOST_AUTO_TEST_CASE(PDNS_Interaction)
{
    namespace test = mjolnir::test;
    mjolnir::LoggerManager::set_default_logger("test_pdns_interaction.log");

    using real_type         = double;
    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;
    using topology_type     = mjolnir::Topology;

    using interaction_type       = mjolnir::ProteinDNANonSpecificInteraction<traits_type>;
    using potential_type         = typename interaction_type::potential_type;
    using parameter_list_type    = typename interaction_type::parameter_list_type;
    using contact_parameter_type = typename parameter_list_type::contact_parameter_type;
    using dna_index_type         = typename parameter_list_type::dna_index_type;
    using ignore_molecule_type   = typename parameter_list_type::ignore_molecule_type;
    using ignore_group_type      = typename parameter_list_type::ignore_group_type;
    using partition_type         = mjolnir::NaivePairCalculation<traits_type, potential_type>;

    //  0 o              |
    //     \             |
    //    1 o === o 3    |
    //     /       \     |
    //  2 o         o 4  |

    constexpr auto pi  = mjolnir::math::constants<real_type>::pi();

    const real_type r0     = 5.0;
    const real_type theta0 = pi * 0.5;
    const real_type phi0   = pi * 2.0 / 3.0;

    const real_type sigma  = 1.0;
    const real_type delta  = pi / 18.0;

    const contact_parameter_type p_pro{1, 0, 2, -1.2, r0, theta0, phi0, r0+5.0, (r0+5.0)*(r0+5.0)};
    const dna_index_type p_dna{3, 4};

    parameter_list_type parameter_list(
        sigma, delta, 5.0, {p_pro}, {p_dna}, {/*exclusion*/},
        ignore_molecule_type("Nothing"), ignore_group_type({}));

    interaction_type interaction(potential_type{},
        parameter_list_type(parameter_list),
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()));

    topology_type topol(5);
    topol.add_connection(0, 1, "bond");
    topol.add_connection(1, 2, "bond");
    topol.add_connection(3, 4, "bond");
    topol.construct_molecules();

    system_type sys(5, boundary_type{});
    test::clear_everything(sys);

    interaction.initialize(sys, topol);

    std::mt19937 rng(123456789);

    for(const auto r     : {r0 - 0.2, r0 + 0.5})
    {
    for(const auto theta : {theta0 + delta * 0.2,
                            theta0 + delta * 0.7,
                            theta0 + delta * 1.2,
                            theta0 - delta * 1.2,
                            theta0 - delta * 0.7,
                            theta0 - delta * 0.2})
    {
    for(const auto phi   : {phi0 + delta * 0.2,
                            phi0 + delta * 0.7,
                            phi0 + delta * 1.2,
                            phi0 - delta * 1.2,
                            phi0 - delta * 0.7,
                            phi0 - delta * 0.2})
    {
        BOOST_TEST_MESSAGE("=======================================================");
        BOOST_TEST_MESSAGE("r = " << r << ", theta = " << theta << ", phi = " << phi);
        BOOST_TEST_MESSAGE("theta0 = " << theta0 << ", delta = " << delta);

        test::clear_everything(sys);

        //  2 o              |
        //     \             |
        //  ^ 1 o === o 3    |
        //     /       \     |
        //  0 o         o 4  |

        sys.position(1) = mjolnir::math::make_coordinate<coord_type>(0.0, 0.0, 0.0);
        sys.position(0) = mjolnir::math::make_coordinate<coord_type>(
                -std::cos(theta) - 1.0, -std::sin(theta), 0.0);
        sys.position(2) = mjolnir::math::make_coordinate<coord_type>(
                 std::cos(theta) - 1.0,  std::sin(theta), 0.0);

        sys.position(3) = mjolnir::math::make_coordinate<coord_type>(
                r, 0.0, 0.0);
        sys.position(4) = mjolnir::math::make_coordinate<coord_type>(
                r + std::cos(pi-phi), -std::sin(pi-phi), 0.0);

        test::apply_random_rotation(sys, rng);

        // check configuration is okay
        {
            // XXX here it does not use boundary condition.
            const auto v13 = sys.position(3) - sys.position(1);
            BOOST_TEST_REQUIRE(mjolnir::math::length(v13) == r,
                               boost::test_tools::tolerance(1e-3));

            const auto v02 = sys.position(2) - sys.position(0);
            const auto dot_1302 = mjolnir::math::dot_product(v02, v13);
            const auto cos_1302 = mjolnir::math::clamp<real_type>(dot_1302 *
                mjolnir::math::rlength(v13) * mjolnir::math::rlength(v02), -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos_1302) == theta,
                               boost::test_tools::tolerance(1e-3));

            const auto v43 = sys.position(3) - sys.position(4);
            const auto dot_134 = mjolnir::math::dot_product(v13, v43);
            const auto cos_134 = mjolnir::math::clamp<real_type>(dot_134 *
                mjolnir::math::rlength(v13) * mjolnir::math::rlength(v43), -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos_134) == phi,
                               boost::test_tools::tolerance(1e-3));
        }

        test::apply_random_perturbation(sys, rng, 0.01);

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, interaction, tol, dr);
        test::check_virial(sys, interaction, tol);
        test::check_force_and_virial(sys, interaction, tol);
        test::check_force_and_energy(sys, interaction, tol);

    } // theta2
    } // theta1
    } // rbp
}
