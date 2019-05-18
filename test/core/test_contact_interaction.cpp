#define BOOST_TEST_MODULE "test_contact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/ContactInteraction.hpp>
#include <mjolnir/potential/local/GoContactPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(Contact_calc_force)
{
    mjolnir::LoggerManager::set_default_logger("test_Contact_calc_force");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::GoContactPotential<real_type>;
    using interaction_type = mjolnir::ContactInteraction<traits_type, potential_type>;

    constexpr real_type tol = 1e-8;

    auto normalize = [](const coord_type& v){return v / mjolnir::math::length(v);};

    const real_type k(100.);
    const real_type native(2.0);

    potential_type   potential(k, native);
    interaction_type interaction("none", {{ {{0,1}}, potential}});

    system_type sys(2, boundary_type{});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;

    sys.at(0).position = coord_type(0,0,0);
    sys.at(1).position = coord_type(0,0,0);
    sys.at(0).velocity = coord_type(0,0,0);
    sys.at(1).velocity = coord_type(0,0,0);
    sys.at(0).force    = coord_type(0,0,0);
    sys.at(1).force    = coord_type(0,0,0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";

    const real_type dr = 1e-3;
    real_type dist = 1e0;
    sys[1].position[0] = dist;

    interaction.initialize(sys);

    for(int i = 0; i < 2000; ++i)
    {
        sys[0].position = coord_type(0,0,0);
        sys[1].position = coord_type(0,0,0);
        sys[0].force    = coord_type(0,0,0);
        sys[1].force    = coord_type(0,0,0);
        sys[1].position[0] = dist;

        interaction.update_margin(dr, sys);

        const real_type deriv = potential.derivative(dist);
        const real_type coef  = std::abs(deriv);

        interaction.calc_force(sys);

        const real_type force_strength1 = mjolnir::math::length(sys.force(0));
        const real_type force_strength2 = mjolnir::math::length(sys.force(1));
        BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
        BOOST_TEST(coef == force_strength2, boost::test_tools::tolerance(tol));

        // direction
        if(dist < native) // repulsive
        {
            const real_type dir1 = mjolnir::math::dot_product(
                normalize(sys[0].force), normalize(sys[0].position - sys[1].position));
            const real_type dir2 = mjolnir::math::dot_product(
                normalize(sys[1].force), normalize(sys[1].position - sys[0].position));

            BOOST_TEST(dir1 == 1.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dir2 == 1.0, boost::test_tools::tolerance(tol));
        }
        else
        {
            const real_type dir1 = mjolnir::math::dot_product(
                normalize(sys[0].force), normalize(sys[1].position - sys[0].position));
            const real_type dir2 = mjolnir::math::dot_product(
                normalize(sys[1].force), normalize(sys[0].position - sys[1].position));

            BOOST_TEST(dir1 == 1e0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dir2 == 1e0, boost::test_tools::tolerance(tol));
        }
        BOOST_TEST(mjolnir::math::length(sys[0].force + sys[1].force) == 0.0,
                   boost::test_tools::tolerance(tol));

        dist += dr;
    }
}
