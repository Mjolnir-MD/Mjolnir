#define BOOST_TEST_MODULE "test_bond_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/BondAngleInteraction.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(BondAngleInteraction_force)
{
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = traits_type::real_type;
    using coord_type      = traits_type::coordinate_type;
    using boundary_type   = traits_type::boundary_type;
    using system_type     = mjolnir::System<traits_type>;
    using harmonic_type   = mjolnir::HarmonicPotential<real_type>;
    using bond_angle_type = mjolnir::BondAngleInteraction<traits_type, harmonic_type>;

    constexpr real_type tol = 1e-7;

    auto normalize = [](const coord_type& v){return v / mjolnir::math::length(v);};

    const real_type k(1e0);
    const real_type native(mjolnir::math::constants<real_type>::pi() * 2.0 / 3.0); // 120 degree

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    const coord_type pos1(1., 0., 0.);
    const coord_type pos2(0., 0., 0.);
    system_type sys(3, boundary_type{});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(2).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;
    sys.at(2).rmass = 1.0;

    sys.at(0).position = pos1;
    sys.at(1).position = pos2;
    sys.at(2).position = coord_type(0,0,0);
    sys.at(0).velocity = coord_type(0,0,0);
    sys.at(1).velocity = coord_type(0,0,0);
    sys.at(2).velocity = coord_type(0,0,0);
    sys.at(0).force    = coord_type(0,0,0);
    sys.at(1).force    = coord_type(0,0,0);
    sys.at(2).force    = coord_type(0,0,0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(2).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";
    sys.at(2).group = "NONE";

    const int N = 1800;
    const real_type dtheta = mjolnir::math::constants<real_type>::pi()  / N;
    for(int i = 1; i < N; ++i)
    {
        BOOST_TEST(mjolnir::math::length(sys[0].position - pos1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[1].position - pos2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[0].velocity)        == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[1].velocity)        == 0.0, boost::test_tools::tolerance(tol));
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
        BOOST_TEST(mjolnir::math::length(sys[0].position - sys[1].position) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[2].position - sys[1].position) == 1.0, boost::test_tools::tolerance(tol));

        const real_type force_strength1 = mjolnir::math::length(sys[0].force);
        const real_type force_strength3 = mjolnir::math::length(sys[2].force);
        if(i == 1200) // most stable point
        {
            BOOST_TEST(coef == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == 0.0, boost::test_tools::tolerance(tol));
        }
        else
        {
            BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == force_strength3, boost::test_tools::tolerance(tol));
        }

        // force applied to center particle is equal to sum of others
        const coord_type sum = sys[0].force + sys[1].force + sys[2].force;
        BOOST_TEST(mjolnir::math::length(sum) == 0.0, boost::test_tools::tolerance(tol));

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_TEST(std::abs(force_strength1) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(std::abs(force_strength3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i < 1200) // narrow
        {
            // perpendicular to radius vector
            const real_type dot1 = mjolnir::math::dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = mjolnir::math::dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., -1., 0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[0].force) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta + mjolnir::math::constants<real_type>::pi() / 2.),
                    sin(theta + mjolnir::math::constants<real_type>::pi() / 2.),
                    0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[2].force) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i > 1200) // extensive
        {
            const real_type dot1 = mjolnir::math::dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = mjolnir::math::dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., 1., 0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[0].force) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta - mjolnir::math::constants<real_type>::pi() / 2.),
                    sin(theta - mjolnir::math::constants<real_type>::pi() / 2.),
                    0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[2].force) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }

        // perpendicular to z axis
        BOOST_TEST(sys[0].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[1].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[2].force[2] == 0.0, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(BondAngleInteraction_numerical_diff)
{
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = traits_type::real_type;
    using coord_type      = traits_type::coordinate_type;
    using boundary_type   = traits_type::boundary_type;
    using system_type     = mjolnir::System<traits_type>;
    using harmonic_type   = mjolnir::HarmonicPotential<real_type>;
    using bond_angle_type = mjolnir::BondAngleInteraction<traits_type, harmonic_type>;

    const real_type k(10.0);
    const real_type native(mjolnir::math::constants<real_type>::pi() / 3.0);

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    for(int i = 0; i < 1000; ++i)
    {
        system_type sys(3, boundary_type{});

        sys.mass (0)= 1.0;
        sys.mass (1)= 1.0;
        sys.mass (2)= 1.0;
        sys.rmass(0) = 1.0;
        sys.rmass(1) = 1.0;
        sys.rmass(2) = 1.0;

        sys.position(0) = coord_type(1.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt));
        sys.position(1) = coord_type(0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt));
        sys.position(2) = coord_type(1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt));
        sys.velocity(0) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(1) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(2) = coord_type(0.0, 0.0, 0.0);
        sys.force   (0) = coord_type(0.0, 0.0, 0.0);
        sys.force   (1) = coord_type(0.0, 0.0, 0.0);
        sys.force   (2) = coord_type(0.0, 0.0, 0.0);

        sys.name (0) = "X";
        sys.name (1) = "X";
        sys.name (2) = "X";
        sys.group(0) = "NONE";
        sys.group(1) = "NONE";
        sys.group(2) = "NONE";

        const auto init = sys;

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;
        for(std::size_t idx=0; idx<3; ++idx)
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

                BOOST_TEST(-dE == dr * mjolnir::math::X(sys.force(idx)),
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

                BOOST_TEST(-dE == dr * mjolnir::math::Y(sys.force(idx)),
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

                BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }

        // -----------------------------------------------------------------
        // check virial
        using matrix33_type = typename traits_type::matrix33_type;

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

BOOST_AUTO_TEST_CASE(BondAngleInteraction_calc_force_and_energy)
{
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = traits_type::real_type;
    using coord_type      = traits_type::coordinate_type;
    using boundary_type   = traits_type::boundary_type;
    using system_type     = mjolnir::System<traits_type>;
    using harmonic_type   = mjolnir::HarmonicPotential<real_type>;
    using bond_angle_type = mjolnir::BondAngleInteraction<traits_type, harmonic_type>;

    const real_type k(10.0);
    const real_type native(mjolnir::math::constants<real_type>::pi() / 3.0);

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    for(int i = 0; i < 1000; ++i)
    {
        system_type sys(3, boundary_type{});

        sys.mass (0)= 1.0;
        sys.mass (1)= 1.0;
        sys.mass (2)= 1.0;
        sys.rmass(0) = 1.0;
        sys.rmass(1) = 1.0;
        sys.rmass(2) = 1.0;

        sys.position(0) = coord_type(1.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt));
        sys.position(1) = coord_type(0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt));
        sys.position(2) = coord_type(1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt));
        sys.velocity(0) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(1) = coord_type(0.0, 0.0, 0.0);
        sys.velocity(2) = coord_type(0.0, 0.0, 0.0);
        sys.force   (0) = coord_type(0.0, 0.0, 0.0);
        sys.force   (1) = coord_type(0.0, 0.0, 0.0);
        sys.force   (2) = coord_type(0.0, 0.0, 0.0);

        sys.name (0) = "X";
        sys.name (1) = "X";
        sys.name (2) = "X";
        sys.group(0) = "NONE";
        sys.group(1) = "NONE";
        sys.group(2) = "NONE";

        constexpr real_type tol = 1e-4;
        auto ref_sys = sys;

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
