#define BOOST_TEST_MODULE "test_omp_dihedral_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/math/math.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/PeriodicGaussianPotential.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/DihedralAngleInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_DihedralAngle_calc_force)
{
    namespace test = mjolnir::test;

    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_dihedral_angle_interaction.log");

    using seq_traits_type  = mjolnir::SimulatorTraits      <double, mjolnir::UnlimitedBoundary>;
    using omp_traits_type  = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type        = typename omp_traits_type::real_type;
    using coordinate_type  = typename omp_traits_type::coordinate_type;
    using boundary_type    = typename omp_traits_type::boundary_type;

    using potential_type   = mjolnir::PeriodicGaussianPotential<real_type>;

    using omp_system_type      = mjolnir::System<omp_traits_type>;
    using omp_interaction_type = mjolnir::DihedralAngleInteraction<omp_traits_type, potential_type>;

    using seq_system_type      = mjolnir::System<seq_traits_type>;
    using seq_interaction_type = mjolnir::DihedralAngleInteraction<seq_traits_type, potential_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    std::mt19937 rng(123456789);

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        omp_system_type omp_sys(10, boundary_type{});
        test::clear_everything(omp_sys);

        omp_sys.position(0) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  0.0,  0.0);
        omp_sys.position(1) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  5.0,  0.0);
        omp_sys.position(2) = mjolnir::math::make_coordinate<coordinate_type>( 5.0,  5.0,  0.0);
        omp_sys.position(3) = mjolnir::math::make_coordinate<coordinate_type>( 5.0,  5.0,  5.0);
        omp_sys.position(4) = mjolnir::math::make_coordinate<coordinate_type>( 5.0, 10.0,  5.0);
        omp_sys.position(5) = mjolnir::math::make_coordinate<coordinate_type>(10.0, 10.0,  5.0);
        omp_sys.position(6) = mjolnir::math::make_coordinate<coordinate_type>(10.0, 10.0, 10.0);
        omp_sys.position(7) = mjolnir::math::make_coordinate<coordinate_type>(10.0, 15.0, 10.0);
        omp_sys.position(8) = mjolnir::math::make_coordinate<coordinate_type>(15.0, 15.0, 10.0);
        omp_sys.position(9) = mjolnir::math::make_coordinate<coordinate_type>(15.0, 15.0, 15.0);

        test::apply_random_perturbation(omp_sys, rng, 0.1);

        // init sequential one with the same coordinates
        seq_system_type seq_sys(10, boundary_type{});
        for(std::size_t i=0; i<seq_sys.size(); ++i)
        {
            seq_sys.mass(i)     = omp_sys.mass(i);
            seq_sys.position(i) = omp_sys.position(i);
            seq_sys.velocity(i) = omp_sys.velocity(i);
            seq_sys.force(i)    = omp_sys.force(i);
            seq_sys.name(i)     = omp_sys.name(i);
            seq_sys.group(i)    = omp_sys.group(i);
        }

        potential_type   potential(/* k = */-100.0, /*sigma = */ 3.1415 * 0.5,
                                   /* native angle = */ -3.1415 * 0.5);

        omp_interaction_type omp_interaction(/*topol = */"none", {
                {std::array<std::size_t, 4>{{0, 1, 2, 3}}, potential},
                {std::array<std::size_t, 4>{{1, 2, 3, 4}}, potential},
                {std::array<std::size_t, 4>{{2, 3, 4, 5}}, potential},
                {std::array<std::size_t, 4>{{3, 4, 5, 6}}, potential},
                {std::array<std::size_t, 4>{{4, 5, 6, 7}}, potential},
                {std::array<std::size_t, 4>{{5, 6, 7, 8}}, potential},
                {std::array<std::size_t, 4>{{6, 7, 8, 9}}, potential}
            });

        seq_interaction_type seq_interaction("none", {
                {std::array<std::size_t, 4>{{0, 1, 2, 3}}, potential},
                {std::array<std::size_t, 4>{{1, 2, 3, 4}}, potential},
                {std::array<std::size_t, 4>{{2, 3, 4, 5}}, potential},
                {std::array<std::size_t, 4>{{3, 4, 5, 6}}, potential},
                {std::array<std::size_t, 4>{{4, 5, 6, 7}}, potential},
                {std::array<std::size_t, 4>{{5, 6, 7, 8}}, potential},
                {std::array<std::size_t, 4>{{6, 7, 8, 9}}, potential}
            });

        omp_interaction.initialize(omp_sys);
        seq_interaction.initialize(seq_sys);

        test::check_force_consistency              (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_and_energy_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_and_virial_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_energy_virial_consistency(omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_energy_consistency             (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
    }
}
