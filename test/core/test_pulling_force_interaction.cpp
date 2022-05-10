#define BOOST_TEST_MODULE "test_pulling_force_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/external/PullingForceInteraction.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(PullingForce)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_pulling_force_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coordinate_type  = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using interaction_type = mjolnir::PullingForceInteraction<traits_type>;

    interaction_type interaction(
            std::vector<std::pair<std::size_t, coordinate_type>>{
                {0, coordinate_type( 0.1,  0.2,  0.3)},
                {2, coordinate_type(-0.3, -0.4, -0.5)}
            });

    system_type sys(4, boundary_type{});

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)  = 1.0;
        sys.rmass(i) = 1.0;

        sys.position(i) = coordinate_type( 0.0,  0.0,  0.0);
        sys.velocity(i) = coordinate_type( 0.0,  0.0,  0.0);
        sys.force   (i) = coordinate_type( 0.0,  0.0,  0.0);

        sys.name (i) = "X";
        sys.group(i) = "NONE";
    }

    interaction.calc_force(sys);

    const real_type tol = 1e-6;

    BOOST_TEST( 0.1 == mjolnir::math::X(sys.force(0)), boost::test_tools::tolerance(tol));
    BOOST_TEST( 0.2 == mjolnir::math::Y(sys.force(0)), boost::test_tools::tolerance(tol));
    BOOST_TEST( 0.3 == mjolnir::math::Z(sys.force(0)), boost::test_tools::tolerance(tol));

    BOOST_TEST( 0.0 == mjolnir::math::X(sys.force(1)), boost::test_tools::tolerance(tol));
    BOOST_TEST( 0.0 == mjolnir::math::Y(sys.force(1)), boost::test_tools::tolerance(tol));
    BOOST_TEST( 0.0 == mjolnir::math::Z(sys.force(1)), boost::test_tools::tolerance(tol));

    BOOST_TEST(-0.3 == mjolnir::math::X(sys.force(2)), boost::test_tools::tolerance(tol));
    BOOST_TEST(-0.4 == mjolnir::math::Y(sys.force(2)), boost::test_tools::tolerance(tol));
    BOOST_TEST(-0.5 == mjolnir::math::Z(sys.force(2)), boost::test_tools::tolerance(tol));

    BOOST_TEST( 0.0 == mjolnir::math::X(sys.force(3)), boost::test_tools::tolerance(tol));
    BOOST_TEST( 0.0 == mjolnir::math::Y(sys.force(3)), boost::test_tools::tolerance(tol));
    BOOST_TEST( 0.0 == mjolnir::math::Z(sys.force(3)), boost::test_tools::tolerance(tol));
}
