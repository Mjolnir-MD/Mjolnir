#define BOOST_TEST_MODULE "test_bond_angle_interaction"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/constants.hpp>
#include <mjolnir/potential/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(BondAngleInteraction_force)
{
    typedef mjolnir::SimulatorTraitsBase<double, mjolnir::UnlimitedBoundary> traits;
    constexpr static traits::real_type tolerance = 1e-7;

    typedef traits::real_type real_type;
    typedef traits::coordinate_type            coord_type;
    typedef traits::boundary_type              boundary_type;
    typedef mjolnir::Particle<coord_type>      particle_type;
    typedef mjolnir::System<traits>            system_type;
    typedef mjolnir::HarmonicPotential<traits> harmonic_type;
    typedef mjolnir::BondAngleInteraction<traits, harmonic_type> bond_angle_type;

    auto normalize = [](const coord_type& v){return v / mjolnir::length(v);};

    const real_type k(1e0);
    const real_type native(mjolnir::constants<real_type>::pi * 2.0 / 3.0); // 120 degree

    harmonic_type potential{k, native};
    bond_angle_type interaction({{ {{0,1,2}}, potential}});

    const coord_type pos1(1., 0., 0.);
    const coord_type pos2(0., 0., 0.);
    std::vector<particle_type> ps{
        {1., pos1,              coord_type(0,0,0), coord_type(0,0,0)},
        {1., pos2,              coord_type(0,0,0), coord_type(0,0,0)},
        {1., coord_type(0,0,0), coord_type(0,0,0), coord_type(0,0,0)}
    };
    system_type sys(std::move(ps), boundary_type{});

    const std::size_t N = 1800;
    const real_type dtheta = mjolnir::constants<real_type>::pi  / N;
    for(int i = 1; i < N; ++i)
    {
        BOOST_CHECK_SMALL(length(sys[0].position - pos1), tolerance);
        BOOST_CHECK_SMALL(length(sys[1].position - pos2), tolerance);
        BOOST_CHECK_SMALL(length(sys[0].velocity), tolerance);
        BOOST_CHECK_SMALL(length(sys[1].velocity), tolerance);
        sys[0].force = coord_type(0,0,0);
        sys[1].force = coord_type(0,0,0);
        sys[2].force = coord_type(0,0,0);

        const real_type theta = i * dtheta;
        const coord_type pos3(std::cos(theta), std::sin(theta), 0e0);
        sys[2].position = pos3;

        const real_type deriv = potential.derivative(theta);
        const real_type coef = std::abs(deriv);

        interaction.calc_force(sys);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV/dtheta.
        BOOST_CHECK_CLOSE_FRACTION(length(sys[0].position - sys[1].position), 1e0, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(length(sys[2].position - sys[1].position), 1e0, tolerance);

        const real_type force_strength1 = length(sys[0].force);
        const real_type force_strength3 = length(sys[2].force);
        if(i == 1200) // most stable point
        {
            BOOST_CHECK_SMALL(coef, tolerance);
            BOOST_CHECK_SMALL(coef, tolerance);
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength1, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(coef, force_strength3, tolerance);
        }

        // force applied to center particle is equal to sum of others
        const coord_type sum = sys[0].force + sys[1].force + sys[2].force;
        BOOST_CHECK_SMALL(length(sum), tolerance);

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_CHECK_SMALL(std::abs(force_strength1), tolerance);
            BOOST_CHECK_SMALL(std::abs(force_strength3), tolerance);
        }
        else if(i < 1200) // narrow
        {
            // perpendicular to radius vector
            const real_type dot1 = dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_CHECK_SMALL(dot1, tolerance);
            BOOST_CHECK_SMALL(dot2, tolerance);

            const coord_type f1(0., -1., 0.);
            BOOST_CHECK_SMALL(length(normalize(sys[0].force) - f1), tolerance);

            const coord_type f3(
                    cos(theta + mjolnir::constants<real_type>::pi / 2.),
                    sin(theta + mjolnir::constants<real_type>::pi / 2.),
                    0.);
            BOOST_CHECK_SMALL(length(normalize(sys[2].force) - f3), tolerance);
        }
        else if(i > 1200) // extensive
        {
            const real_type dot1 = dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_CHECK_SMALL(dot1, tolerance);
            BOOST_CHECK_SMALL(dot2, tolerance);

            const coord_type f1(0., 1., 0.);
            BOOST_CHECK_SMALL(length(normalize(sys[0].force) - f1), tolerance);

            const coord_type f3(
                    cos(theta - mjolnir::constants<real_type>::pi / 2.),
                    sin(theta - mjolnir::constants<real_type>::pi / 2.),
                    0.);
            BOOST_CHECK_SMALL(length(normalize(sys[2].force) - f3), tolerance);
        }

        // perpendicular to z axis
        BOOST_CHECK_SMALL(sys[0].force[2], tolerance);
        BOOST_CHECK_SMALL(sys[1].force[2], tolerance);
        BOOST_CHECK_SMALL(sys[2].force[2], tolerance);
    }
}

// BOOST_AUTO_TEST_CASE(BondAngleInteraction_energy)
// {
//     const real_type k(1e0);
//     const real_type native(M_PI * 2.0 / 3.0);
//     LocalPotentialBaseSptr potential =
//         std::make_shared<HarmonicPotential>(k, native);
//
//     InteractionBaseSptr<3> angl(new BondAngleInteraction(potential));
//
//     const coord_type pos1(1e0, 0e0, 0e0);
//     const coord_type pos2(0e0, 0e0, 0e0);
//     const real_type dtheta = M_PI / 1800.0;
//     for(int i = 0; i < 1000; ++i)
//     {
//         const real_type theta = i * dtheta;
//         const coord_type pos3(std::cos(theta), std::sin(theta), 0e0);
//
//         const real_type pot = potential->calc_potential(theta);
//
//         const real_type ene = angl->calc_energy(
//                 std::array<coord_type, 3>{{pos1, pos2, pos3}});
//
//         if(i == 1200) // native
//             BOOST_CHECK_SMALL(ene, tolerance);
//         else 
//             BOOST_CHECK_CLOSE_FRACTION(pot, ene, tolerance);
//     }
// }
//

