#define BOOST_TEST_MODULE "test_position_restraint_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/traits.hpp>
#include <mjolnir/interaction/external/PositionRestraintInteraction.hpp>
#include <mjolnir/potential/external/HarmonicRestraintPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(PositionRestraint_Harmonic)
{
    typedef mjolnir::test::traits<double> traits;
    constexpr static traits::real_type tol = 1e-8;

    using real_type        = traits::real_type;
    using coordinate_type  = traits::coordinate_type;
    using boundary_type    = traits::boundary_type;
    using system_type      = mjolnir::System<traits>;
    using potential_type   = mjolnir::HarmonicRestraintPotential<real_type>;
    using interaction_type = mjolnir::PositionRestraintInteraction<traits, potential_type>;

    auto normalize = [](const coordinate_type& v){return v / mjolnir::math::length(v);};

    potential_type   potential(std::vector<std::pair<real_type, real_type>>{
            {/*k =*/1.0, /*r0 =*/ 0.0},
            {/*k =*/1.0, /*r0 =*/10.0},
        });
    interaction_type interaction(
        std::vector<std::pair<std::size_t, coordinate_type>>{
            {0, coordinate_type{0.0,0.0,0.0}},
            {1, coordinate_type{0.0,0.0,0.0}}
        }, potential);

    system_type sys(2, boundary_type{});

    sys.at(0).mass  = 1.0;
    sys.at(1).mass  = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;

    sys.at(0).position = coordinate_type( 0.0,  0.0,  0.0);
    sys.at(1).position = coordinate_type(10.0,  0.0,  0.0);
    sys.at(0).velocity = coordinate_type( 0.0,  0.0,  0.0);
    sys.at(1).velocity = coordinate_type( 0.0,  0.0,  0.0);
    sys.at(0).force    = coordinate_type( 0.0,  0.0,  0.0);
    sys.at(1).force    = coordinate_type( 0.0,  0.0,  0.0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";

    std::mt19937 rng(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);
    std::normal_distribution<real_type> gauss(0.0, 1.0);

    for(int i = 0; i < 10000; ++i)
    {
        sys.position(0) = coordinate_type(uni(rng), uni(rng), uni(rng));
        sys.force(0)    = coordinate_type(0, 0, 0);

        const auto dist = mjolnir::math::length(sys.position(0));
        const auto deriv = potential.derivative(0, dist);
        const auto coef  = std::abs(deriv);

        interaction.calc_force(sys);

        const auto force_strength = mjolnir::math::length(sys.force(0));

        BOOST_TEST(coef == force_strength, boost::test_tools::tolerance(tol));

        // forces always attract the particle to the origin.
        const auto direction = mjolnir::math::dot_product(
            normalize(sys.force(0)), normalize(-sys.position(0)));

        BOOST_TEST(direction == 1e0, boost::test_tools::tolerance(tol));
    }

    for(int i = 0; i < 10000; ++i)
    {
        sys.position(1) = coordinate_type(uni(rng), uni(rng), uni(rng));
        sys.force(1)    = coordinate_type(0, 0, 0);

        const auto offset = 10.0 * normalize(
                coordinate_type(gauss(rng), gauss(rng), gauss(rng)));
        sys.position(1) += offset; // dist ~ 10.0 +- random vector

        const auto dist  = mjolnir::math::length(sys.position(1));
        const auto deriv = potential.derivative(1, dist);
        const auto coef  = std::abs(deriv);

        interaction.calc_force(sys);
        const auto force_strength = mjolnir::math::length(sys.force(1));

        BOOST_TEST(coef == force_strength, boost::test_tools::tolerance(tol));

        if(dist < 10.0) // repulsive
        {
            const auto direction = mjolnir::math::dot_product(
                normalize(sys.force(1)), normalize(sys.position(1)));
            BOOST_TEST(direction == 1e0, boost::test_tools::tolerance(tol));
        }
        else // attractive
        {
            const auto direction = mjolnir::math::dot_product(
                normalize(sys.force(1)), normalize(-sys.position(1)));
            BOOST_TEST(direction == 1e0, boost::test_tools::tolerance(tol));
        }
    }
}