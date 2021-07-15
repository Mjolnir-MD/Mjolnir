#define BOOST_TEST_MODULE "test_bond_length_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/BondLengthInteraction.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(BondLength_numerical_difference)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_bondlength_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    const real_type k(100.0);
    const real_type native(std::sqrt(3.0));

    potential_type   potential(k, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    std::mt19937 mt(123456789);

    for(std::size_t i = 0; i < 1000; ++i)
    {
        system_type sys(2, boundary_type{});
        test::clear_everything(sys);

        sys.position(0) = coord_type( 0.0, 0.0, 0.0);
        sys.position(1) = coord_type( 1.0, 1.0, 1.0);

        test::apply_random_rotation(sys, mt);
        test::apply_random_perturbation(sys, mt, 0.01);

        constexpr real_type tol = 1e-6;
        constexpr real_type dr  = 1e-6;

        test::check_force(sys, interaction, tol, dr);
        test::check_virial(sys, interaction, tol);
        test::check_force_and_energy(sys, interaction, tol);
    }
}


BOOST_AUTO_TEST_CASE(BondLength_calc_force)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_bondlength_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type    = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    constexpr real_type tol = 1e-8;

    auto normalize = [](const coord_type& v){return v / mjolnir::math::length(v);};

    const real_type k(100.);
    const real_type native(2.0);

    potential_type    potential(k, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    system_type sys(2, boundary_type{});
    test::clear_everything(sys);

    const real_type dr = 1e-3;
    real_type dist = 1e0;
    for(int i = 0; i < 2000; ++i)
    {
        sys.position(0) = coord_type(0,0,0);
        sys.position(1) = coord_type(0,0,0);
        sys.force(0)    = coord_type(0,0,0);
        sys.force(1)    = coord_type(0,0,0);
        sys.position(1)[0] = dist;

        const real_type deriv = potential.derivative(dist);
        const real_type coef  = std::abs(deriv);

        interaction.calc_force(sys);

        const real_type force_strength1 = mjolnir::math::length(sys.force(0));
        const real_type force_strength2 = mjolnir::math::length(sys.force(1));


        // direction
        if(i == 1000) // most stable point
        {
            BOOST_TEST(force_strength1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(force_strength2 == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i < 1000) // repulsive
        {
            BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == force_strength2, boost::test_tools::tolerance(tol));

            const real_type dir1 = mjolnir::math::dot_product(
                normalize(sys.force(0)), normalize(sys.position(0) - sys.position(1)));
            const real_type dir2 = mjolnir::math::dot_product(
                normalize(sys.force(1)), normalize(sys.position(1) - sys.position(0)));

            BOOST_TEST(dir1 == 1.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dir2 == 1.0, boost::test_tools::tolerance(tol));
        }
        else if(i > 1000) // attractive
        {
            BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == force_strength2, boost::test_tools::tolerance(tol));

            const real_type dir1 = mjolnir::math::dot_product(
                normalize(sys.force(0)), normalize(sys.position(1) - sys.position(0)));
            const real_type dir2 = mjolnir::math::dot_product(
                normalize(sys.force(1)), normalize(sys.position(0) - sys.position(1)));

            BOOST_TEST(dir1 == 1e0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dir2 == 1e0, boost::test_tools::tolerance(tol));
        }
        BOOST_TEST(mjolnir::math::length(sys.force(0) + sys.force(1)) == 0.0,
                   boost::test_tools::tolerance(tol));

        dist += dr;
    }
}
