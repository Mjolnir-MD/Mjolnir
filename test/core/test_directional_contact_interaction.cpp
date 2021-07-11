#define BOOST_TEST_MODULE "test_directional_contact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/DirectionalContactInteraction.hpp>
#include <mjolnir/forcefield/local/CosinePotential.hpp>
#include <mjolnir/forcefield/local/GoContactPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/math/Matrix.hpp>

#include <tuple>
#include <random>

BOOST_AUTO_TEST_CASE(DirectionalContactInteraction_numerical_diff)
{
    namespace test = mjolnir::test;
    mjolnir::LoggerManager::set_default_logger("test_directional_contact_interaction.log");
    using traits_type              = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type                = traits_type::real_type;
    using coord_type               = traits_type::coordinate_type;
    using boundary_type            = traits_type::boundary_type;
    using system_type              = mjolnir::System<traits_type>;
    using angle_potential_type     = mjolnir::CosinePotential<real_type>;
    using contact_potential_type   = mjolnir::GoContactPotential<real_type>;
    using directional_contact_type = mjolnir::DirectionalContactInteraction<
        traits_type, angle_potential_type, angle_potential_type, contact_potential_type>;

    constexpr real_type pi = mjolnir::math::constants<real_type>::pi();
    const     real_type k_angle(1e0);
    const     real_type k_contact(1e0);
    const     real_type angle1_native(pi);
    const     real_type angle2_native(pi);
    const     real_type contact_native(2.0);
    const auto angle1_potential    = angle_potential_type{k_angle, 1, angle1_native};
    const auto angle2_potential    = angle_potential_type{k_angle, 1, angle2_native};
    const auto contact_potential   = contact_potential_type{k_contact, contact_native};
    directional_contact_type interaction("none",
        {{std::make_tuple(std::array<std::size_t, 4>{{0,1,2,3}},
        angle1_potential, angle2_potential, contact_potential)
        }});

    system_type sys(4, boundary_type{});
    sys.mass(0) = 1.0;
    sys.mass(1) = 1.0;
    sys.mass(2) = 1.0;
    sys.mass(3) = 1.0;
    sys.position(0) = mjolnir::math::make_coordinate<coord_type>(0.0, 0.0, 0.0);
    sys.position(1) = mjolnir::math::make_coordinate<coord_type>(0.0, 0.0, 0.0);
    sys.position(2) = mjolnir::math::make_coordinate<coord_type>(0.0, 0.0, 0.0);
    sys.position(3) = mjolnir::math::make_coordinate<coord_type>(0.0, 0.0, 0.0);
    interaction.initialize(sys);

    const int angle_step_num           = 108;
    const std::size_t contact_step_num = 10;
    const real_type dtheta             = pi / real_type(angle_step_num);
    const real_type contact_cutoff     = contact_potential.cutoff();
    const real_type test_contact_range = 1.1 * contact_cutoff;
    const real_type dr                 = test_contact_range / real_type(contact_step_num);
    const real_type tol                = 1e-3;
    const real_type dx                 = 1e-6;

    std::mt19937 mt(123456789);

    for(std::size_t i = 0; i < angle_step_num; ++i)
    {
        for(std::size_t j = 0; j < angle_step_num; ++j)
        {
            for(std::size_t k = 1; k <= contact_step_num; ++k)
            {
                const real_type theta1 = i * dtheta;
                const real_type theta2 = j * dtheta;
                if(std::sin(theta1) < 1e-2 || std::sin(theta2) < 1e-2)
                {
                    // In this case, numerical error become big and difficult to
                    // implement appropriate test. So skip.
                    continue;
                }

                const real_type r = k * dr;
                BOOST_TEST_MESSAGE("theta1 = " << theta1 << ", theta2 = " << pi + theta2
                                   << ", r = " << r);

                sys.position(0) = coord_type(std::cos(theta1), std::sin(theta1), 0e0);
                sys.position(1) = coord_type(0e0, 0e0, 0e0);
                sys.position(2) = coord_type(r, 0e0, 0e0);
                sys.position(3) = coord_type(r - std::cos(theta2), std::sin(theta2), 0e0);
                sys.velocity(0) = coord_type(0.0, 0.0, 0.0);
                sys.velocity(1) = coord_type(0.0, 0.0, 0.0);
                sys.velocity(2) = coord_type(0.0, 0.0, 0.0);
                sys.velocity(3) = coord_type(0.0, 0.0, 0.0);
                sys.force   (0) = coord_type(0.0, 0.0, 0.0);
                sys.force   (1) = coord_type(0.0, 0.0, 0.0);
                sys.force   (2) = coord_type(0.0, 0.0, 0.0);
                sys.force   (3) = coord_type(0.0, 0.0, 0.0);

                test::apply_random_rotation(sys, mt);
                test::apply_random_perturbation(sys, mt, 0.01);

                interaction.update(sys);

                test::check_force(sys, interaction, tol, dx);
                test::check_virial(sys, interaction, tol);
                test::check_force_and_energy(sys, interaction, tol);

                sys.force(0) = coord_type(0.0, 0.0, 0.0);
                sys.force(1) = coord_type(0.0, 0.0, 0.0);
                sys.force(2) = coord_type(0.0, 0.0, 0.0);
                sys.force(3) = coord_type(0.0, 0.0, 0.0);

                const auto total_f = sys.force(0) + sys.force(1) + sys.force(2) + sys.force(3);
                BOOST_TEST_REQUIRE(mjolnir::math::X(total_f) == 0.0, boost::test_tools::tolerance(tol));
                BOOST_TEST_REQUIRE(mjolnir::math::Y(total_f) == 0.0, boost::test_tools::tolerance(tol));
                BOOST_TEST_REQUIRE(mjolnir::math::Z(total_f) == 0.0, boost::test_tools::tolerance(tol));
            }
        }
    }
}
