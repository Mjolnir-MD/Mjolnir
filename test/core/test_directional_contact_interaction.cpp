#define BOOST_TEST_MODULE "test_directional_contact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/DirectionalContactInteraction.hpp>
#include <mjolnir/potential/local/CosinePotential.hpp>
#include <mjolnir/potential/local/GoContactPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/math/Matrix.hpp>

#include <tuple>
#include <random>

BOOST_AUTO_TEST_CASE(DirectionalContactInteraction_numerical_diff)
{
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

  const real_type pi = mjolnir::math::constants<real_type>::pi();
  const real_type k_angle(1e0);
  const real_type k_contact(1e0);
  const real_type angle1_native(pi);
  const real_type angle2_native(pi);
  const real_type contact_native(2.0);
  angle_potential_type angle1_potential    = angle_potential_type{k_angle, 1, angle1_native};
  angle_potential_type angle2_potential    = angle_potential_type{k_angle, 1, angle2_native};
  contact_potential_type contact_potential = contact_potential_type{k_contact, contact_native};
  directional_contact_type interaction("none",
      {{std::make_tuple(std::array<std::size_t, 4>{0,1,2,3},
      angle1_potential, angle2_potential, contact_potential)
      }});

  system_type sys(4, boundary_type{});
  sys.mass(0) = 1.0;
  sys.mass(1) = 1.0;
  sys.mass(2) = 1.0;
  sys.mass(3) = 1.0;
  sys.at(0).position = {0.0, 0.0, 0.0};
  sys.at(1).position = {0.0, 0.0, 0.0};
  sys.at(2).position = {0.0, 0.0, 0.0};
  sys.at(3).position = {0.0, 0.0, 0.0};
  interaction.initialize(sys);

  const int angle_step_num           = 108;
  const std::size_t contact_step_num = 10;
  const real_type dtheta             = pi / real_type(angle_step_num);
  const real_type contact_cutoff     = contact_potential.cutoff();
  const real_type test_contact_range = 1.1 * contact_cutoff;
  const real_type dr                 = test_contact_range / real_type(contact_step_num);
  const real_type tol                = 1e-3;
  const real_type dx                 = 1e-6;
  const real_type ignore_tolerance   = 1e-12;
  std::mt19937 mt(123456789);
  std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

  for(std::size_t i = 0; i < angle_step_num; ++i)
    {
      for(std::size_t j = 0; j < angle_step_num; ++j)
        {
          for(std::size_t k = 1; k <= contact_step_num; ++k)
            {
              real_type theta1 = i * dtheta;
              real_type theta2 = j * dtheta;
              if(std::sin(theta1) < 1e-2 || std::sin(theta2) < 1e-2)
              {
                // In this case, numerical error become big and difficult to
                // implement appropriate test. So skip.
                continue;
              }

              real_type r = k * dr;
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


              // rotate randomly system
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

              const auto init = sys;
              interaction.update(sys);

              for(std::size_t idx=0; idx<4; ++idx)
                {
                  BOOST_TEST_MESSAGE("idx = " << idx);
                  {
                    sys = init;
                    interaction.calc_force(sys);
                    const auto total_f = sys.force(0) + sys.force(1) + sys.force(2) + sys.force(3);
                    BOOST_TEST_REQUIRE(mjolnir::math::X(total_f) == 0e0,
                                       boost::test_tools::tolerance(tol));
                    BOOST_TEST_REQUIRE(mjolnir::math::Y(total_f) == 0e0,
                                       boost::test_tools::tolerance(tol));
                    BOOST_TEST_REQUIRE(mjolnir::math::Z(total_f) == 0e0,
                                       boost::test_tools::tolerance(tol));
                  }
                  { // test for X coordinate
                    // reset positions
                    sys = init;

                    // calc F(x)
                    interaction.calc_force(sys);

                    // calcu U(x-dx)
                    mjolnir::math::X(sys.position(idx)) -= dx;
                    const auto E0 = interaction.calc_energy(sys);

                    // calc U(x+dx)
                    mjolnir::math::X(sys.position(idx)) += 2.0 * dx;
                    const auto E1 = interaction.calc_energy(sys);

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;
                    if(std::abs(dE) < ignore_tolerance)
                    {
                      BOOST_TEST(dE == 0.0, boost::test_tools::tolerance(tol));
                    }
                    else
                    {
                      BOOST_TEST(-dE == dx * mjolnir::math::X(sys.force(idx)),
                                 boost::test_tools::tolerance(tol));
                    }
                  }
                  { // test for Y coordinate
                    // reset positions
                    sys = init;

                    // calc F(x)
                    interaction.calc_force(sys);

                    // calc U(x-dx)
                    mjolnir::math::Y(sys.position(idx)) -= dx;
                    const auto E0 = interaction.calc_energy(sys);

                    // calc U(x+dx)
                    mjolnir::math::Y(sys.position(idx)) += 2.0 * dx;
                    const auto E1 = interaction.calc_energy(sys);

                    //central difference
                    const auto dE = (E1 - E0) * 0.5;
                    if(std::abs(dE) < ignore_tolerance)
                    {
                      BOOST_TEST(dE == 0.0, boost::test_tools::tolerance(tol));
                    }
                    else
                    {
                      BOOST_TEST(-dE == dx * mjolnir::math::Y(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                    }
                  }
                  { // test for Z coordinate
                    // reset positions
                    sys = init;

                    // calc F(x)
                    interaction.calc_force(sys);

                    // calc U(x-dx)
                    mjolnir::math::Z(sys.position(idx)) -= dx;
                    const auto E0 = interaction.calc_energy(sys);

                    // calc U(x+dx)
                    mjolnir::math::Z(sys.position(idx)) += 2.0 * dx;
                    const auto E1 = interaction.calc_energy(sys);

                    //central difference
                    const auto dE = (E1 - E0) * 0.5;

                    if(std::abs(dE) < ignore_tolerance)
                    {
                      BOOST_TEST(dE == 0.0, boost::test_tools::tolerance(tol));
                    }
                    else
                    {
                      BOOST_TEST(-dE == dx * mjolnir::math::Z(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                    }
                  }
                }
            }
        }
    }
}
