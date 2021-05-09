#define BOOST_TEST_MODULE "test_bond_length_gocontact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/BondLengthGoContactInteraction.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(BondLengthGoContactInteraction)
{
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::GoContactPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    constexpr real_type tol = 1e-7;

    const coord_type pos1(1.0, 0.0, 0.0);
    const coord_type pos2(0.0, 0.0, 0.0);
    system_type sys(2, boundary_type{});

    sys.at(0).mass     = 1.0;
    sys.at(1).mass     = 1.0;
    sys.at(0).rmass    = 1.0;
    sys.at(1).rmass    = 1.0;
    sys.at(0).position = pos1;
    sys.at(1).position = pos2;
    sys.at(0).velocity = coord_type(0,0,0);
    sys.at(1).velocity = coord_type(0,0,0);
    sys.at(0).force    = coord_type(0,0,0);
    sys.at(1).force    = coord_type(0,0,0);
    sys.at(0).name     = "X";
    sys.at(1).name     = "X";
    sys.at(0).group    = "NONE";
    sys.at(1).group    = "NONE";

    const real_type k(1.0);
    const real_type native(1.0);

    potential_type   potential(k, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    const int N = 10000;
    const real_type dx = 0.001;
    for(int i = 1; i < N; ++i)
    {
        const auto dist = 0.6 + i * dx;
        sys.position(0) = coord_type(0.0,  0.0, 0.0);
        sys.position(1) = coord_type(dist, 0.0, 0.0);
        sys.force(0)    = coord_type(0.0, 0.0, 0.0);
        sys.force(1)    = coord_type(0.0, 0.0, 0.0);

        const real_type deriv = potential.derivative(dist);

        interaction.calc_force(sys);

        const real_type force_strength1 = mjolnir::math::length(sys.force(0));
        const real_type force_strength2 = mjolnir::math::length(sys.force(1));

        BOOST_TEST(force_strength1 == std::abs(deriv), boost::test_tools::tolerance(tol));
        BOOST_TEST(force_strength2 == std::abs(deriv), boost::test_tools::tolerance(tol));

        const auto E_pot = potential.potential(dist);
        const auto E_int = interaction.calc_energy(sys);

        BOOST_TEST(E_pot == E_int, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(BondLengthGoContact_numerical_difference)
{
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::GoContactPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    const real_type k(1.0);
    const real_type native(std::sqrt(3.0));

    potential_type   potential(k, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    for(std::size_t i = 0; i < 1000; ++i)
    {
        system_type sys(2, boundary_type{});

        sys.at(0).mass  = 1.0;
        sys.at(1).mass  = 1.0;
        sys.at(0).rmass = 1.0;
        sys.at(1).rmass = 1.0;

        sys.at(0).position = coord_type( 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt));
        sys.at(1).position = coord_type( 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt));
        sys.at(0).velocity = coord_type( 0.0, 0.0, 0.0);
        sys.at(1).velocity = coord_type( 0.0, 0.0, 0.0);
        sys.at(0).force    = coord_type( 0.0, 0.0, 0.0);
        sys.at(1).force    = coord_type( 0.0, 0.0, 0.0);

        sys.at(0).name  = "X";
        sys.at(1).name  = "X";
        sys.at(0).group = "TEST";
        sys.at(1).group = "TEST";

        const auto init = sys;

        constexpr real_type tol = 2e-4;
        constexpr real_type dr  = 1e-5;
        for(std::size_t idx=0; idx<2; ++idx)
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

BOOST_AUTO_TEST_CASE(BondLengthGoContact_calc_force_and_energy)
{
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::GoContactPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    const real_type k(1.0);
    const real_type native(std::sqrt(3.0));

    potential_type   potential(k, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    for(std::size_t i = 0; i < 1000; ++i)
    {
        system_type sys(2, boundary_type{});

        sys.at(0).mass  = 1.0;
        sys.at(1).mass  = 1.0;
        sys.at(0).rmass = 1.0;
        sys.at(1).rmass = 1.0;

        sys.at(0).position = coord_type( 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt));
        sys.at(1).position = coord_type( 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt));
        sys.at(0).velocity = coord_type( 0.0, 0.0, 0.0);
        sys.at(1).velocity = coord_type( 0.0, 0.0, 0.0);
        sys.at(0).force    = coord_type( 0.0, 0.0, 0.0);
        sys.at(1).force    = coord_type( 0.0, 0.0, 0.0);

        sys.at(0).name  = "X";
        sys.at(1).name  = "X";
        sys.at(0).group = "TEST";
        sys.at(1).group = "TEST";

        auto ref_sys = sys;
        constexpr real_type tol = 2e-4;

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
