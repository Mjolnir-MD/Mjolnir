#define BOOST_TEST_MODULE "test_PWMcos_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/PWMcos/PWMcosInteraction.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>

BOOST_AUTO_TEST_CASE(PWMcos_Interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_pdns_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;
    using topology_type     = mjolnir::Topology;

    using real_type              = double;
    using potential_type         = mjolnir::PWMcosPotential<traits_type>;
    using base_kind              = typename potential_type::base_kind;
    using contact_parameter_type = typename potential_type::contact_parameter_type;
    using dna_parameter_type     = typename potential_type::dna_parameter_type;
    using interaction_type       = mjolnir::PWMcosInteraction<traits_type>;
    using ignore_molecule_type   = typename potential_type::ignore_molecule_type;
    using ignore_group_type      = typename potential_type::ignore_group_type;
    using partition_type         = mjolnir::NaivePairCalculation<traits_type, potential_type>;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    // protein     DNA
    //
    //  0 o       3 o --o  |
    //     \     4         |
    //    1 o === o -- o 6 |
    //     /               |
    //  2 o       5 o--o   |

    constexpr auto pi  = mjolnir::math::constants<real_type>::pi();

    const real_type r0       = 5.0;
    const real_type theta1_0 = pi * 0.4;
    const real_type theta2_0 = pi * 0.5;
    const real_type theta3_0 = pi * 0.6;

    const real_type sigma  = 1.0;
    const real_type delta  = pi / 18.0;

    const contact_parameter_type p_pro{1, 0, 2, -1.2, 0.0, r0, theta1_0, theta2_0, theta3_0, r0 + 5*sigma, {{1.0, 1.0, 1.0, 1.0}} };
    const dna_parameter_type     p_dna{base_kind::A, 4, 6, 3, 5};

    potential_type potential(
        sigma, delta, 1.0, 0.0, 5.0, {p_pro}, {p_dna}, {/*exclusion*/},
        ignore_molecule_type("Nothing"), ignore_group_type({}));

    interaction_type interaction(potential_type(potential),
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()));

    system_type sys(7, boundary_type{});
    topology_type topol(7);

    topol.add_connection(0, 1, "bond");
    topol.add_connection(1, 2, "bond");
    topol.add_connection(3, 4, "bond");
    topol.add_connection(4, 5, "bond");
    topol.add_connection(4, 6, "bond");
    topol.construct_molecules();

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)  = 1.0;
        sys.rmass(i) = 1.0;
        sys.position(i) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(i) = coord_type(0.0, 0.0, 0.0);
        sys.force(i)    = coord_type(0.0, 0.0, 0.0);
        sys.name(i)  = "A";
        sys.group(i) = "A";
    }
    potential  .initialize(sys, topol);
    interaction.initialize(sys, topol);

    for(const auto r     : {r0 - 0.2, r0 + 0.5})
    {
    for(const auto theta1 : {theta1_0 + delta * 0.2,
                             theta1_0 + delta * 0.7,
                             theta1_0 + delta * 1.2,
                             theta1_0 - delta * 1.2,
                             theta1_0 - delta * 0.7,
                             theta1_0 - delta * 0.2})
    {
    for(const auto theta2 : {theta2_0 + delta * 0.2,
                             theta2_0 + delta * 0.7,
                             theta2_0 + delta * 1.2,
                             theta2_0 - delta * 1.2,
                             theta2_0 - delta * 0.7,
                             theta2_0 - delta * 0.2})
    {
    for(const auto theta3 : {theta3_0 + delta * 0.2,
                             theta3_0 + delta * 0.7,
                             theta3_0 + delta * 1.2,
                             theta3_0 - delta * 1.2,
                             theta3_0 - delta * 0.7,
                             theta3_0 - delta * 0.2})
    {
        // r
        sys.position(1) = mjolnir::math::make_coordinate<coord_type>(r, 0, 0);
        sys.position(4) = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        // theta1
        sys.position(6) = mjolnir::math::make_coordinate<coord_type>(std::cos(theta1), std::sin(theta1), 0);
        // theta2
        sys.position(3) = mjolnir::math::make_coordinate<coord_type>(std::cos(theta2+pi), std::sin(theta2+pi), 0);
        sys.position(5) = mjolnir::math::make_coordinate<coord_type>(std::cos(theta2),    std::sin(theta2),    0);
        // theta3
        sys.position(0) = mjolnir::math::make_coordinate<coord_type>(r + std::cos(theta3),    std::sin(theta3),    0);
        sys.position(2) = mjolnir::math::make_coordinate<coord_type>(r + std::cos(theta3+pi), std::sin(theta3+pi), 0);

        // rotate in random direction
        {
            using matrix33_type = typename traits_type::matrix33_type;

            const auto rot_x = uni(mt) * pi;
            const auto rot_y = uni(mt) * pi;
            const auto rot_z = uni(mt) * pi;

            matrix33_type rotm_x(1.0,             0.0,              0.0,
                                 0.0, std::cos(rot_x), -std::sin(rot_x),
                                 0.0, std::sin(rot_x),  std::cos(rot_x));
            matrix33_type rotm_y( std::cos(rot_y), 0.0,  std::sin(rot_y),
                                              0.0, 1.0,              0.0,
                                 -std::sin(rot_y), 0.0,  std::cos(rot_y));
            matrix33_type rotm_z(std::cos(rot_z), -std::sin(rot_z), 0.0,
                                 std::sin(rot_z),  std::cos(rot_z), 0.0,
                                             0.0,              0.0, 1.0);

            const matrix33_type rotm = rotm_x * rotm_y * rotm_z;
            for(std::size_t idx=0; idx<sys.size(); ++idx)
            {
                sys.position(idx) = rotm * sys.position(idx);
            }
        }
        // check configuration is okay
        {
            // XXX here it does not use boundary condition.
            const auto v41 = sys.position(1) - sys.position(4); // B -> C
            BOOST_TEST_REQUIRE(mjolnir::math::length(v41) == r,
                               boost::test_tools::tolerance(1e-3));

            const auto v46 = sys.position(6) - sys.position(4); // B -> S
            const auto v35 = sys.position(5) - sys.position(3); // B-1 -> B+1
            const auto v20 = sys.position(0) - sys.position(2); // C+1 -> C-1

            const auto dot1 = mjolnir::math::dot_product(v41, v46) *
                              mjolnir::math::rlength(v41) *
                              mjolnir::math::rlength(v46);
            const auto cos1 = mjolnir::math::clamp<real_type>(dot1, -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos1) == theta1,
                               boost::test_tools::tolerance(1e-3));

            const auto dot2 = mjolnir::math::dot_product(v41, v35) *
                              mjolnir::math::rlength(v41) *
                              mjolnir::math::rlength(v35);
            const auto cos2 = mjolnir::math::clamp<real_type>(dot2, -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos2) == theta2,
                               boost::test_tools::tolerance(1e-3));

            const auto dot3 = mjolnir::math::dot_product(v41, v20) *
                              mjolnir::math::rlength(v41) *
                              mjolnir::math::rlength(v20);
            const auto cos3 = mjolnir::math::clamp<real_type>(dot3, -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos3) == theta3,
                               boost::test_tools::tolerance(1e-3));
        }

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.position(idx) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.force(idx)     = coord_type(0.0, 0.0, 0.0);
        }
        const system_type init = sys;
        interaction.initialize(init, topol);

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
                BOOST_TEST_MESSAGE("force.x = " << mjolnir::math::X(sys.force(idx)) <<
                                   ", -dE/dr = " << -dE/dr);
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
                BOOST_TEST_MESSAGE("force.y = " << mjolnir::math::Y(sys.force(idx)) <<
                                   ", -dE/dr = " << -dE/dr);
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
                BOOST_TEST_MESSAGE("force.z = " << mjolnir::math::Z(sys.force(idx)) <<
                                   ", -dE/dr = " << -dE/dr);
            }
        }
    } // theta3
    } // theta2
    } // theta1
    } // rbp
}

BOOST_AUTO_TEST_CASE(PWMcos_Interaction_force_and_energy)
{
    mjolnir::LoggerManager::set_default_logger("test_pdns_interaction.log");

    using traits_type       = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type         = traits_type::real_type;
    using coord_type        = traits_type::coordinate_type;
    using boundary_type     = traits_type::boundary_type;
    using system_type       = mjolnir::System<traits_type>;
    using topology_type     = mjolnir::Topology;

    using real_type              = double;
    using potential_type         = mjolnir::PWMcosPotential<traits_type>;
    using base_kind              = typename potential_type::base_kind;
    using contact_parameter_type = typename potential_type::contact_parameter_type;
    using dna_parameter_type     = typename potential_type::dna_parameter_type;
    using interaction_type       = mjolnir::PWMcosInteraction<traits_type>;
    using ignore_molecule_type   = typename potential_type::ignore_molecule_type;
    using ignore_group_type      = typename potential_type::ignore_group_type;
    using partition_type         = mjolnir::NaivePairCalculation<traits_type, potential_type>;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    // protein     DNA
    //
    //  0 o       3 o --o  |
    //     \     4         |
    //    1 o === o -- o 6 |
    //     /               |
    //  2 o       5 o--o   |

    constexpr auto pi  = mjolnir::math::constants<real_type>::pi();

    const real_type r0       = 5.0;
    const real_type theta1_0 = pi * 0.4;
    const real_type theta2_0 = pi * 0.5;
    const real_type theta3_0 = pi * 0.6;

    const real_type sigma  = 1.0;
    const real_type delta  = pi / 18.0;

    const contact_parameter_type p_pro{1, 0, 2, -1.2, 0.0, r0, theta1_0, theta2_0, theta3_0, r0 + 5*sigma, {{1.0, 1.0, 1.0, 1.0}} };
    const dna_parameter_type     p_dna{base_kind::A, 4, 6, 3, 5};

    potential_type potential(
        sigma, delta, 1.0, 0.0, 5.0, {p_pro}, {p_dna}, {/*exclusion*/},
        ignore_molecule_type("Nothing"), ignore_group_type({}));

    interaction_type interaction(potential_type(potential),
        mjolnir::SpatialPartition<traits_type, potential_type>(
            mjolnir::make_unique<partition_type>()));

    system_type sys(7, boundary_type{});
    topology_type topol(7);

    topol.add_connection(0, 1, "bond");
    topol.add_connection(1, 2, "bond");
    topol.add_connection(3, 4, "bond");
    topol.add_connection(4, 5, "bond");
    topol.add_connection(4, 6, "bond");
    topol.construct_molecules();

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)  = 1.0;
        sys.rmass(i) = 1.0;
        sys.position(i) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(i) = coord_type(0.0, 0.0, 0.0);
        sys.force(i)    = coord_type(0.0, 0.0, 0.0);
        sys.name(i)  = "A";
        sys.group(i) = "A";
    }
    potential  .initialize(sys, topol);
    interaction.initialize(sys, topol);

    for(const auto r     : {r0 - 0.2, r0 + 0.5})
    {
    for(const auto theta1 : {theta1_0 + delta * 0.2,
                             theta1_0 + delta * 0.7,
                             theta1_0 + delta * 1.2,
                             theta1_0 - delta * 1.2,
                             theta1_0 - delta * 0.7,
                             theta1_0 - delta * 0.2})
    {
    for(const auto theta2 : {theta2_0 + delta * 0.2,
                             theta2_0 + delta * 0.7,
                             theta2_0 + delta * 1.2,
                             theta2_0 - delta * 1.2,
                             theta2_0 - delta * 0.7,
                             theta2_0 - delta * 0.2})
    {
    for(const auto theta3 : {theta3_0 + delta * 0.2,
                             theta3_0 + delta * 0.7,
                             theta3_0 + delta * 1.2,
                             theta3_0 - delta * 1.2,
                             theta3_0 - delta * 0.7,
                             theta3_0 - delta * 0.2})
    {
        // r
        sys.position(1) = mjolnir::math::make_coordinate<coord_type>(r, 0, 0);
        sys.position(4) = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        // theta1
        sys.position(6) = mjolnir::math::make_coordinate<coord_type>(std::cos(theta1), std::sin(theta1), 0);
        // theta2
        sys.position(3) = mjolnir::math::make_coordinate<coord_type>(std::cos(theta2+pi), std::sin(theta2+pi), 0);
        sys.position(5) = mjolnir::math::make_coordinate<coord_type>(std::cos(theta2),    std::sin(theta2),    0);
        // theta3
        sys.position(0) = mjolnir::math::make_coordinate<coord_type>(r + std::cos(theta3),    std::sin(theta3),    0);
        sys.position(2) = mjolnir::math::make_coordinate<coord_type>(r + std::cos(theta3+pi), std::sin(theta3+pi), 0);

        // rotate in random direction
        {
            using matrix33_type = typename traits_type::matrix33_type;

            const auto rot_x = uni(mt) * pi;
            const auto rot_y = uni(mt) * pi;
            const auto rot_z = uni(mt) * pi;

            matrix33_type rotm_x(1.0,             0.0,              0.0,
                                 0.0, std::cos(rot_x), -std::sin(rot_x),
                                 0.0, std::sin(rot_x),  std::cos(rot_x));
            matrix33_type rotm_y( std::cos(rot_y), 0.0,  std::sin(rot_y),
                                              0.0, 1.0,              0.0,
                                 -std::sin(rot_y), 0.0,  std::cos(rot_y));
            matrix33_type rotm_z(std::cos(rot_z), -std::sin(rot_z), 0.0,
                                 std::sin(rot_z),  std::cos(rot_z), 0.0,
                                             0.0,              0.0, 1.0);

            const matrix33_type rotm = rotm_x * rotm_y * rotm_z;
            for(std::size_t idx=0; idx<sys.size(); ++idx)
            {
                sys.position(idx) = rotm * sys.position(idx);
            }
        }
        // check configuration is okay
        {
            // XXX here it does not use boundary condition.
            const auto v41 = sys.position(1) - sys.position(4); // B -> C
            BOOST_TEST_REQUIRE(mjolnir::math::length(v41) == r,
                               boost::test_tools::tolerance(1e-3));

            const auto v46 = sys.position(6) - sys.position(4); // B -> S
            const auto v35 = sys.position(5) - sys.position(3); // B-1 -> B+1
            const auto v20 = sys.position(0) - sys.position(2); // C+1 -> C-1

            const auto dot1 = mjolnir::math::dot_product(v41, v46) *
                              mjolnir::math::rlength(v41) *
                              mjolnir::math::rlength(v46);
            const auto cos1 = mjolnir::math::clamp<real_type>(dot1, -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos1) == theta1,
                               boost::test_tools::tolerance(1e-3));

            const auto dot2 = mjolnir::math::dot_product(v41, v35) *
                              mjolnir::math::rlength(v41) *
                              mjolnir::math::rlength(v35);
            const auto cos2 = mjolnir::math::clamp<real_type>(dot2, -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos2) == theta2,
                               boost::test_tools::tolerance(1e-3));

            const auto dot3 = mjolnir::math::dot_product(v41, v20) *
                              mjolnir::math::rlength(v41) *
                              mjolnir::math::rlength(v20);
            const auto cos3 = mjolnir::math::clamp<real_type>(dot3, -1, 1);
            BOOST_TEST_REQUIRE(std::acos(cos3) == theta3,
                               boost::test_tools::tolerance(1e-3));
        }

        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.position(idx) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.force(idx)     = coord_type(0.0, 0.0, 0.0);
        }
        system_type ref_sys = sys;
        interaction.initialize(sys, topol);

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

    } // theta3
    } // theta2
    } // theta1
    } // rbp
}
