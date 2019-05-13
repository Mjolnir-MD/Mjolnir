#define BOOST_TEST_MODULE "test_omp_bond_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/math/math.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/BondAngleInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_BondAngle)
{
    constexpr double tol = 1e-8;
    mjolnir::LoggerManager::set_default_logger("test_omp_bond_angle_interaction.log");

    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = typename traits_type::real_type;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::BondAngleInteraction<traits_type, potential_type>;
    using rng_type         = mjolnir::RandomNumberGenerator<traits_type>;

    using sequencial_system_type      = mjolnir::System<
        mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>>;
    using sequencial_interaction_type = mjolnir::BondAngleInteraction<
        mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>, potential_type>;

    const int max_number_of_threads = omp_get_max_threads();
    BOOST_TEST_WARN(max_number_of_threads > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    for(int num_thread=1; num_thread<=max_number_of_threads; ++num_thread)
    {
        omp_set_num_threads(num_thread);
        BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

        potential_type   potential(/* k = */100.0, /* native angle = */ 3.1415 * 0.5);
        interaction_type interaction(/*topol = */"none", {
                {std::array<std::size_t, 3>{{0, 1, 2}}, potential},
                {std::array<std::size_t, 3>{{1, 2, 3}}, potential},
                {std::array<std::size_t, 3>{{2, 3, 4}}, potential},
                {std::array<std::size_t, 3>{{3, 4, 5}}, potential},
                {std::array<std::size_t, 3>{{4, 5, 6}}, potential},
                {std::array<std::size_t, 3>{{5, 6, 7}}, potential},
                {std::array<std::size_t, 3>{{6, 7, 8}}, potential},
                {std::array<std::size_t, 3>{{7, 8, 9}}, potential}
            });

        sequencial_interaction_type seq_interaction("none", {
                {std::array<std::size_t, 3>{{0, 1, 2}}, potential},
                {std::array<std::size_t, 3>{{1, 2, 3}}, potential},
                {std::array<std::size_t, 3>{{2, 3, 4}}, potential},
                {std::array<std::size_t, 3>{{3, 4, 5}}, potential},
                {std::array<std::size_t, 3>{{4, 5, 6}}, potential},
                {std::array<std::size_t, 3>{{5, 6, 7}}, potential},
                {std::array<std::size_t, 3>{{6, 7, 8}}, potential},
                {std::array<std::size_t, 3>{{7, 8, 9}}, potential}
            });

        rng_type    rng(123456789);
        system_type sys(10, boundary_type{});
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.mass(i)     = 1.0;
            sys.velocity(i) = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.force(i)    = mjolnir::math::make_coordinate<coordinate_type>(0, 0, 0);
            sys.name(i)     = "X";
            sys.group(i)    = "TEST";
        }

        sys.position(0) = mjolnir::math::make_coordinate<coordinate_type>( 0.0,  0.0, 0.0);
        sys.position(1) = mjolnir::math::make_coordinate<coordinate_type>( 5.0,  0.0, 0.0);
        sys.position(2) = mjolnir::math::make_coordinate<coordinate_type>( 5.0,  5.0, 0.0);
        sys.position(3) = mjolnir::math::make_coordinate<coordinate_type>(10.0,  5.0, 0.0);
        sys.position(4) = mjolnir::math::make_coordinate<coordinate_type>(10.0, 10.0, 0.0);
        sys.position(5) = mjolnir::math::make_coordinate<coordinate_type>(15.0, 10.0, 0.0);
        sys.position(6) = mjolnir::math::make_coordinate<coordinate_type>(15.0, 15.0, 0.0);
        sys.position(7) = mjolnir::math::make_coordinate<coordinate_type>(20.0, 15.0, 0.0);
        sys.position(8) = mjolnir::math::make_coordinate<coordinate_type>(20.0, 20.0, 0.0);
        sys.position(9) = mjolnir::math::make_coordinate<coordinate_type>(25.0, 20.0, 0.0);

        // add perturbation
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            mjolnir::math::X(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Y(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
            mjolnir::math::Z(sys.position(i)) += rng.uniform_real(-0.1, 0.1);
        }

        // init sequential one with the same coordinates
        sequencial_system_type seq_sys(10, boundary_type{});
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            seq_sys.mass(i)     = sys.mass(i);
            seq_sys.position(i) = sys.position(i);
            seq_sys.velocity(i) = sys.velocity(i);
            seq_sys.force(i)    = sys.force(i);
            seq_sys.name(i)     = sys.name(i);
            seq_sys.group(i)    = sys.group(i);
        }

#pragma omp parallel
        {
            // calculate forces with openmp
            interaction.calc_force(sys);
        }
        sys.merge_forces();

        // calculate forces without openmp
        seq_interaction.calc_force(seq_sys);

        // check the values are the same
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            BOOST_TEST(mjolnir::math::X(seq_sys.force(i)) == mjolnir::math::X(sys.force(i)),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Y(seq_sys.force(i)) == mjolnir::math::Y(sys.force(i)),
                       boost::test_tools::tolerance(tol));
            BOOST_TEST(mjolnir::math::Z(seq_sys.force(i)) == mjolnir::math::Z(sys.force(i)),
                       boost::test_tools::tolerance(tol));
        }

        BOOST_TEST(interaction.calc_energy(sys) == seq_interaction.calc_energy(seq_sys),
                   boost::test_tools::tolerance(tol));
    }
}
