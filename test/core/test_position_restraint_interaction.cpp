#define BOOST_TEST_MODULE "test_position_restraint_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/external/PositionRestraintInteraction.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(PositionRestraint_Harmonic)
{
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coordinate_type  = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::PositionRestraintInteraction<traits_type, potential_type>;

    constexpr real_type tol = 1e-8;

    auto normalize = [](const coordinate_type& v){return v / mjolnir::math::length(v);};

    const potential_type pot1{/*k =*/1.0, /*r0 =*/ 0.0};
    const potential_type pot2{/*k =*/1.0, /*r0 =*/10.0};

    interaction_type interaction(
        std::vector<std::tuple<std::size_t, coordinate_type, potential_type>>{
            {0, coordinate_type{0.0,0.0,0.0}, pot1},
            {1, coordinate_type{0.0,0.0,0.0}, pot2}
        });

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
        const auto deriv = pot1.derivative(dist);
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
        const auto deriv = pot2.derivative(dist);
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
