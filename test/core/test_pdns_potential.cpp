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
    using ignore_group_type    = typename potential_type::ignore_group_type;
    using ignore_molecule_type = typename potential_type::ignore_molecule_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi();

    potential_type potential(
        1.0, pi / 18.0, 5.0, {/* contact = empty */}, {/*dnas = empty*/},
        {/* exclude = empty */}, ignore_molecule_type("Nothing"), ignore_group_type({}));

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
    using ignore_group_type    = typename potential_type::ignore_group_type;
    using ignore_molecule_type = typename potential_type::ignore_molecule_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi();

    potential_type potential(
        1.0, pi / 18.0, 5.0, {/* contact = empty */}, {/*dnas = empty*/},
        {/* exclude = empty */}, ignore_molecule_type("Nothing"), ignore_group_type({}));

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
