#define BOOST_TEST_MODULE "test_PDNS_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/System.hpp>
#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificPotential.hpp>

BOOST_AUTO_TEST_CASE(PDNS_Potential_f)
{
    mjolnir::LoggerManager::set_default_logger("test_pdns_potential.log");

    using real_type = double;
    using potential_type = mjolnir::ProteinDNANonSpecificPotential<real_type>;
    using ignore_molecule_type = typename potential_type::ignore_molecule_type;
    using ignore_group_type    = typename potential_type::ignore_group_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi();

    mjolnir::ProteinDNANonSpecificPotential<real_type> potential(
        1.0, pi / 18.0, 5.0, {/* parameter = empty */}, {/* exclude = empty */},
        ignore_molecule_type("Nothing"), ignore_group_type({}));

    const real_type r_0   = 7.0;
    const real_type r_min = 1.0;
    const real_type r_max = 15.0;
    const real_type dr    = (r_max - r_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type r   = r_min + i * dr;
        const real_type f1  = potential.f(r_0, r + h);
        const real_type f2  = potential.f(r_0, r - h);
        const real_type df_numeric  = (f1 - f2) / (2 * h);
        const real_type df_analytic = potential.f_df(r_0, r).second;

        if((f1 == 0.0 || f2 == 0.0) && !(f1 == 0.0 && f2 == 0.0))
        {
            // the numeric differentiation becomes unstable here.
        }
        else
        {
            BOOST_TEST(df_numeric == df_analytic, boost::test_tools::tolerance(h));
        }

        const real_type f1_ = potential.f_df(r_0, r+h).first;
        const real_type f2_ = potential.f_df(r_0, r-h).first;
        BOOST_TEST(f1 == f1_, boost::test_tools::tolerance(h));
        BOOST_TEST(f2 == f2_, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(PDNS_Potential_g)
{
    mjolnir::LoggerManager::set_default_logger("test_pdns_potential.log");

    using real_type = double;
    using potential_type = mjolnir::ProteinDNANonSpecificPotential<real_type>;
    using ignore_molecule_type = typename potential_type::ignore_molecule_type;
    using ignore_group_type    = typename potential_type::ignore_group_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi();

    mjolnir::ProteinDNANonSpecificPotential<real_type> potential(
        1.0, pi / 18.0, 5.0, {/* parameter = empty */}, {/* exclude = empty */},
        ignore_molecule_type("Nothing"), ignore_group_type({}));

    const real_type theta_0   = pi * 0.5;
    const real_type theta_min = pi * 0.3;
    const real_type theta_max = pi * 0.7;
    const real_type dtheta    = (theta_max - theta_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type theta = theta_min + i * dtheta;
        const real_type g1  = potential.g(theta_0, theta + h);
        const real_type g2  = potential.g(theta_0, theta - h);
        const real_type dg_numeric  = (g1 - g2) / (2 * h);
        const real_type dg_analytic = potential.g_dg(theta_0, theta).second;

        if((g1 == 0.0 || g2 == 0.0) && !(g1 == 0.0 && g2 == 0.0))
        {
            // the numeric differentiation becomes unstable here.
        }
        else
        {
            BOOST_TEST(dg_numeric == dg_analytic, boost::test_tools::tolerance(h));
        }
        const real_type g1_ = potential.g_dg(theta_0, theta+h).first;
        const real_type g2_ = potential.g_dg(theta_0, theta-h).first;
        BOOST_TEST(g1 == g1_, boost::test_tools::tolerance(h));
        BOOST_TEST(g2 == g2_, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(PDNS_Potential_interacting_pair)
{
    mjolnir::LoggerManager::set_default_logger("test_pdns_potential.log");

    using real_type = double;
    using potential_type = mjolnir::ProteinDNANonSpecificPotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;
    using bead_kind      = typename potential_type::bead_kind;
    using ignore_molecule_type = typename potential_type::ignore_molecule_type;
    using ignore_group_type    = typename potential_type::ignore_group_type;

    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi();
    const auto nil = potential_type::invalid();

    const parameter_type p_pro{bead_kind::Protein, nil, 1,   2,   1.0, 2.0, 3.0, 4.0};
    const parameter_type p_dna{bead_kind::DNA,     3,   nil, nil, 0.0, 0.0, 0.0, 0.0};

    mjolnir::ProteinDNANonSpecificPotential<real_type> potential(
        1.0, pi / 18.0, 5.0, {{0, p_pro}, {1, p_pro}, {2, p_dna}, {3, p_dna}},
        {/* exclude */}, ignore_molecule_type("Nothing"), ignore_group_type({})
    );

    using traits_type     = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type   = typename traits_type::boundary_type;
    using coordinate_type = typename traits_type::coordinate_type;

    mjolnir::System<traits_type> sys(4, boundary_type{});
    for(std::size_t i=0; i<4; ++i)
    {
        sys.mass(i)     = 0.0;
        sys.rmass(i)    = 0.0;
        sys.position(i) = mjolnir::math::make_coordinate<coordinate_type>(0,0,0);
        sys.velocity(i) = mjolnir::math::make_coordinate<coordinate_type>(0,0,0);
        sys.force(i)    = mjolnir::math::make_coordinate<coordinate_type>(0,0,0);
        sys.name(i)     = "X";
        sys.group(i)    = "NONE";
    }

    sys.topology().construct_molecules();
    potential.initialize(sys);

    BOOST_TEST(!potential.has_interaction(0, 1));
    BOOST_TEST( potential.has_interaction(0, 2));
    BOOST_TEST( potential.has_interaction(0, 3));
    BOOST_TEST( potential.has_interaction(1, 2));
    BOOST_TEST( potential.has_interaction(1, 3));
    BOOST_TEST(!potential.has_interaction(2, 3));

    const auto pp = potential.prepare_params(0, 2);
    BOOST_TEST(pp.PN == 1);
    BOOST_TEST(pp.PC == 2);
    BOOST_TEST(pp.S3 == 3);

    BOOST_TEST(pp.k      == 1.0);
    BOOST_TEST(pp.r0     == 2.0);
    BOOST_TEST(pp.theta0 == 3.0);
    BOOST_TEST(pp.phi0   == 4.0);
}
