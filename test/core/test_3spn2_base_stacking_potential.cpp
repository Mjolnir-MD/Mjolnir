#define BOOST_TEST_MODULE "test_3SPN2_BaseStacking_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/local/ThreeSPN2BaseStackingPotential.hpp>

BOOST_AUTO_TEST_CASE(f_3SPN2_BaseStacking_double)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = double;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi;

    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    const real_type theta_0   = 0.5  * pi;
    const real_type theta_min = 0.01 * pi;
    const real_type theta_max = 0.99 * pi;
    const real_type dtheta    = (theta_max - theta_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type theta = theta_min + i * dtheta;
        const real_type f1  = potential.f(theta + h, theta_0);
        const real_type f2  = potential.f(theta - h, theta_0);
        const real_type df_numeric  = (f1 - f2) / (2 * h);
        const real_type df_analytic = potential.df(theta, theta_0);

        if((f1 == 0.0 || f2 == 0.0) && !(f1 == 0.0 && f2 == 0.0))
        {
            // the numeric differentiation becomes unstable.
        }
        else
        {
            BOOST_TEST(df_numeric == df_analytic, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(f_3SPN2_BaseStacking_float)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = float;

    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3f;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi;

    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    const real_type theta_0  = 0.5  * pi;
    const real_type theta_min = 0.01 * pi;
    const real_type theta_max = 0.99 * pi;
    const real_type dtheta    = (theta_max - theta_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type theta = theta_min + i * dtheta;
        const real_type f1  = potential.f(theta + h, theta_0);
        const real_type f2  = potential.f(theta - h, theta_0);
        const real_type df_numeric  = (f1 - f2) / (2 * h);
        const real_type df_analytic = potential.df(theta, theta_0);

        if((f1 == 0.0f || f2 == 0.0f) && !(f1 == 0.0f && f2 == 0.0f))
        {
            // the numeric differentiation becomes unstable.
        }
        else
        {
            BOOST_TEST(df_numeric == df_analytic, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(U_rep_3SPN2_BaseStacking_double)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = double;
    using base_stack_kind = mjolnir::parameter_3SPN2::base_stack_kind;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    const real_type x_min = 2.0;
    const real_type x_max = 6.0; // outside of cutoff
    const real_type dx    = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x  = x_min + i * dx;
        const real_type U1 = potential.U_rep(base_stack_kind::AA, x + h);
        const real_type U2 = potential.U_rep(base_stack_kind::AA, x - h);
        const real_type dU_numeric  = (U1 - U2) / (2 * h);
        const real_type dU_analytic = potential.dU_rep(base_stack_kind::AA, x);

        if(U1 != real_type(0.0) && U2 != real_type(0.0))
        {
            BOOST_TEST(dU_numeric == dU_analytic, boost::test_tools::tolerance(h));
        }
        // purely repulsive potential.
        BOOST_TEST(U1 <= U2); // U1 is for x+h, U2 is for x-h
    }
}

BOOST_AUTO_TEST_CASE(U_rep_3SPN2_BaseStacking_float)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = float;
    using base_stack_kind = mjolnir::parameter_3SPN2::base_stack_kind;

    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3f;

    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    const real_type x_min = 2.0;
    const real_type x_max = 6.0; // outside of cutoff
    const real_type dx    = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x  = x_min + i * dx;
        const real_type U1 = potential.U_rep(base_stack_kind::AA, x + h);
        const real_type U2 = potential.U_rep(base_stack_kind::AA, x - h);
        const real_type dU_numeric  = (U1 - U2) / (2 * h);
        const real_type dU_analytic = potential.dU_rep(base_stack_kind::AA, x);

        if(U1 != real_type(0.0) && U2 != real_type(0.0))
        {
            BOOST_TEST(dU_numeric == dU_analytic, boost::test_tools::tolerance(h));
        }
        // purely repulsive potential.
        BOOST_TEST(U1 <= U2); // U1 is for x+h, U2 is for x-h
    }
}

BOOST_AUTO_TEST_CASE(U_attr_3SPN2_BaseStacking_double)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = double;
    using base_stack_kind = mjolnir::parameter_3SPN2::base_stack_kind;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    const real_type x_min =  3.0; // inside of flat region
    const real_type x_max = 12.0;
    const real_type dx    = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x  = x_min + i * dx;
        const real_type U1 = potential.U_attr(base_stack_kind::AA, x + h);
        const real_type U2 = potential.U_attr(base_stack_kind::AA, x - h);
        const real_type dU_numeric  = (U1 - U2) / (2 * h);
        const real_type dU_analytic = potential.dU_attr(base_stack_kind::AA, x);

        if(U1 != potential.epsilon(base_stack_kind::AA) &&
           U2 == potential.epsilon(base_stack_kind::AA))
        {
            BOOST_TEST(dU_numeric == dU_analytic, boost::test_tools::tolerance(h));
        }
        // purely attractive potential
        BOOST_TEST(U2 <= U1);
    }
}

BOOST_AUTO_TEST_CASE(U_attr_3SPN2_BaseStacking_float)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = float;
    using base_stack_kind = mjolnir::parameter_3SPN2::base_stack_kind;

    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3f;

    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    const real_type x_min =  3.0f;
    const real_type x_max = 12.0f;
    const real_type dx    = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x  = x_min + i * dx;
        const real_type U1 = potential.U_attr(base_stack_kind::AA, x + h);
        const real_type U2 = potential.U_attr(base_stack_kind::AA, x - h);
        const real_type dU_numeric  = (U1 - U2) / (2 * h);
        const real_type dU_analytic = potential.dU_attr(base_stack_kind::AA, x);

        if(U1 != potential.epsilon(base_stack_kind::AA) &&
           U2 == potential.epsilon(base_stack_kind::AA))
        {
            BOOST_TEST(dU_numeric == dU_analytic, boost::test_tools::tolerance(h));
        }
        // purely attractive potential
        BOOST_TEST(U2 <= U1);

    }
}

BOOST_AUTO_TEST_CASE(U_dU_attr_3SPN2_BaseStacking_double)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = double;
    using base_stack_kind = mjolnir::parameter_3SPN2::base_stack_kind;

    constexpr std::size_t N   =  100;
    constexpr real_type   h   = 1e-6;
    constexpr real_type x_min =  3.0;
    constexpr real_type x_max = 12.0;
    constexpr real_type dx    = (x_max - x_min) / N;
    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    for(std::size_t i=0; i<N; ++i)
    {
        const auto x    = x_min + i * dx;
        const auto U_dU = potential.U_dU_attr(base_stack_kind::AA, x);
        const auto U    = potential.U_attr(base_stack_kind::AA, x);
        const auto dU   = potential.dU_attr(base_stack_kind::AA, x);

        BOOST_TEST(U_dU.first  ==  U, boost::test_tools::tolerance(h));
        BOOST_TEST(U_dU.second == dU, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(U_dU_attr_3SPN2_BaseStacking_float)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = float;
    using base_stack_kind = mjolnir::parameter_3SPN2::base_stack_kind;

    constexpr std::size_t N   =  100;
    constexpr real_type   h   = 1e-3f;
    constexpr real_type x_min =  3.0f;
    constexpr real_type x_max = 12.0f;
    constexpr real_type dx    = (x_max - x_min) / N;
    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    for(std::size_t i=0; i<N; ++i)
    {
        const auto x    = x_min + i * dx;
        const auto U_dU = potential.U_dU_attr(base_stack_kind::AA, x);
        const auto U    = potential.U_attr(base_stack_kind::AA, x);
        const auto dU   = potential.dU_attr(base_stack_kind::AA, x);

        BOOST_TEST(U_dU.first  ==  U, boost::test_tools::tolerance(h));
        BOOST_TEST(U_dU.second == dU, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(bs_kind_3SPN2_BaseStacking)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_base_stacking_potential.log");
    using real_type = double;
    using base_kind = mjolnir::parameter_3SPN2::base_kind;
    using base_stack_kind = mjolnir::parameter_3SPN2::base_stack_kind;

    mjolnir::ThreeSPN2BaseStackingPotential<real_type> potential;

    BOOST_TEST(base_stack_kind::AA == potential.bs_kind(base_kind::A, base_kind::A));
    BOOST_TEST(base_stack_kind::AT == potential.bs_kind(base_kind::A, base_kind::T));
    BOOST_TEST(base_stack_kind::AG == potential.bs_kind(base_kind::A, base_kind::G));
    BOOST_TEST(base_stack_kind::AC == potential.bs_kind(base_kind::A, base_kind::C));
    BOOST_TEST(base_stack_kind::TA == potential.bs_kind(base_kind::T, base_kind::A));
    BOOST_TEST(base_stack_kind::TT == potential.bs_kind(base_kind::T, base_kind::T));
    BOOST_TEST(base_stack_kind::TG == potential.bs_kind(base_kind::T, base_kind::G));
    BOOST_TEST(base_stack_kind::TC == potential.bs_kind(base_kind::T, base_kind::C));
    BOOST_TEST(base_stack_kind::GA == potential.bs_kind(base_kind::G, base_kind::A));
    BOOST_TEST(base_stack_kind::GT == potential.bs_kind(base_kind::G, base_kind::T));
    BOOST_TEST(base_stack_kind::GG == potential.bs_kind(base_kind::G, base_kind::G));
    BOOST_TEST(base_stack_kind::GC == potential.bs_kind(base_kind::G, base_kind::C));
    BOOST_TEST(base_stack_kind::CA == potential.bs_kind(base_kind::C, base_kind::A));
    BOOST_TEST(base_stack_kind::CT == potential.bs_kind(base_kind::C, base_kind::T));
    BOOST_TEST(base_stack_kind::CG == potential.bs_kind(base_kind::C, base_kind::G));
    BOOST_TEST(base_stack_kind::CC == potential.bs_kind(base_kind::C, base_kind::C));
}
