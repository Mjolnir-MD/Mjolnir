#define BOOST_TEST_MODULE "test_3SPN2BaseBaseInteraction_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteractionPotential.hpp>

BOOST_AUTO_TEST_CASE(f_3SPN2_BaseBaseInteraction_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");

    using real_type = double;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind  ;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi;

    mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type> potential({},
        typename mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>::ignore_group_type({})
        );

    const real_type theta_0   = 0.5  * pi;
    const real_type theta_min = 0.01 * pi;
    const real_type theta_max = 0.99 * pi;
    const real_type dtheta    = (theta_max - theta_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type theta = theta_min + i * dtheta;
        const real_type f1  = potential.f(base_pair_kind::AT, theta + h, theta_0);
        const real_type f2  = potential.f(base_pair_kind::AT, theta - h, theta_0);
        const real_type df_numeric  = (f1 - f2) / (2 * h);
        const real_type df_analytic = potential.df(base_pair_kind::AT, theta, theta_0);

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

BOOST_AUTO_TEST_CASE(f_3SPN2_BaseBaseInteraction_float)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");
    using real_type = float;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind  ;

    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3f;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi;

    mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type> potential({},
        typename mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>::ignore_group_type({})
        );

    const real_type theta_0  = 0.5  * pi;
    const real_type theta_min = 0.01 * pi;
    const real_type theta_max = 0.99 * pi;
    const real_type dtheta    = (theta_max - theta_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type theta = theta_min + i * dtheta;
        const real_type f1  = potential.f(base_pair_kind::AT, theta + h, theta_0);
        const real_type f2  = potential.f(base_pair_kind::AT, theta - h, theta_0);
        const real_type df_numeric  = (f1 - f2) / (2 * h);
        const real_type df_analytic = potential.df(base_pair_kind::AT, theta, theta_0);

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

BOOST_AUTO_TEST_CASE(U_rep_3SPN2_BaseBaseInteraction_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");

    using real_type = double;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind  ;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type> potential({},
        typename mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>::ignore_group_type({})
        );

    const real_type r_min =  3.0;
    const real_type r_max = 12.0;
    const real_type dr    = (r_max - r_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const auto r         = r_min + i * dr;
        const real_type pot1 = potential.U_rep(base_pair_kind::AT, r + h);
        const real_type pot2 = potential.U_rep(base_pair_kind::AT, r - h);
        const real_type df_numeric  = (pot1 - pot2) / (2 * h);
        const real_type df_analytic = potential.dU_rep(base_pair_kind::AT, r);

        if(pot1 != 0.0 && pot2 == 0.0)
        {
            // the numeric differentiation becomes unstable.
        }
        else if(df_analytic == 0.0)
        {
            BOOST_TEST(potential.U_rep(base_pair_kind::AT, r) == 0.0,
                       boost::test_tools::tolerance(h));
        }
        else
        {
            BOOST_TEST(df_numeric == df_analytic, boost::test_tools::tolerance(h));
        }
    }
}
BOOST_AUTO_TEST_CASE(U_attr_3SPN2_BaseBaseInteraction_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");

    using real_type = double;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind  ;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type> potential({},
        typename mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>::ignore_group_type({})
        );

    const real_type r_min =  3.0; // avoiding the long tail around cutoff;
    const real_type r_max = 12.0; // not take numerical error around zero into account
    const real_type dr    = (r_max - r_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const auto r         = r_min + i * dr;
        const real_type pot1 = potential.U_attr(base_pair_kind::AT, r + h);
        const real_type pot2 = potential.U_attr(base_pair_kind::AT, r - h);
        const real_type df_numeric  = (pot1 - pot2) / (2 * h);
        const real_type df_analytic = potential.dU_attr(base_pair_kind::AT, r);

        if(pot1 != 0.0 && pot2 == 0.0)
        {
            // the numeric differentiation becomes unstable.
        }
        else if(df_analytic == 0.0)
        {
            // dU == 0 means that r is in `U == -epsilon` part.
            BOOST_TEST(potential.U_attr(base_pair_kind::AT, r) ==
                       -potential.epsilon(base_pair_kind::AT),
                       boost::test_tools::tolerance(h));
        }
        else
        {
            BOOST_TEST(df_numeric == df_analytic, boost::test_tools::tolerance(1e-5));
        }
    }
}


BOOST_AUTO_TEST_CASE(U_dU_attr_3SPN2_BaseBaseInteraction_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");
    using real_type = double;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind  ;
    using cross_stack_kind = mjolnir::parameter_3SPN2::cross_stack_kind;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-8;

    mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type> potential({},
        typename mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>::ignore_group_type({})
        );
    for(const auto bp : {base_pair_kind::AT, base_pair_kind::TA,
                         base_pair_kind::GC, base_pair_kind::CG})
    {
        const real_type r_min =  3.0;
        const real_type r_max = 18.0;
        const real_type dr    = (r_max - r_min) / N;

        for(std::size_t i=0; i<N; ++i)
        {
            const auto r         = r_min + i * dr;
            const auto U_attr    = potential.U_attr(bp, r);
            const auto dU_attr   = potential.dU_attr(bp, r);
            const auto U_dU_attr = potential.U_dU_attr(bp, r);

            BOOST_TEST(U_dU_attr.first  == U_attr,  boost::test_tools::tolerance(h));
            BOOST_TEST(U_dU_attr.second == dU_attr, boost::test_tools::tolerance(h));
        }
    }
    for(const auto cs : {cross_stack_kind::AA5, cross_stack_kind::AA3,
                         cross_stack_kind::AT5, cross_stack_kind::AT3,
                         cross_stack_kind::AG5, cross_stack_kind::AG3,
                         cross_stack_kind::AC5, cross_stack_kind::AC3,
                         cross_stack_kind::TA5, cross_stack_kind::TA3,
                         cross_stack_kind::TT5, cross_stack_kind::TT3,
                         cross_stack_kind::TG5, cross_stack_kind::TG3,
                         cross_stack_kind::TC5, cross_stack_kind::TC3,
                         cross_stack_kind::GA5, cross_stack_kind::GA3,
                         cross_stack_kind::GT5, cross_stack_kind::GT3,
                         cross_stack_kind::GG5, cross_stack_kind::GG3,
                         cross_stack_kind::GC5, cross_stack_kind::GC3,
                         cross_stack_kind::CA5, cross_stack_kind::CA3,
                         cross_stack_kind::CT5, cross_stack_kind::CT3,
                         cross_stack_kind::CG5, cross_stack_kind::CG3,
                         cross_stack_kind::CC5, cross_stack_kind::CC3})
    {
        const real_type r_min =  3.0;
        const real_type r_max = 18.0;
        const real_type dr    = (r_max - r_min) / N;

        for(std::size_t i=0; i<N; ++i)
        {
            const auto r         = r_min + i * dr;
            const auto U_attr    = potential.U_attr(cs, r);
            const auto dU_attr   = potential.dU_attr(cs, r);
            const auto U_dU_attr = potential.U_dU_attr(cs, r);

            BOOST_TEST(U_dU_attr.first  == U_attr,  boost::test_tools::tolerance(h));
            BOOST_TEST(U_dU_attr.second == dU_attr, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(U_dU_attr_3SPN2_BaseBaseInteraction_float)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");
    using real_type = float;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind;
    using cross_stack_kind = mjolnir::parameter_3SPN2::cross_stack_kind;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-4;

    mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type> potential({},
        typename mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>::ignore_group_type({})
        );

    for(const auto bp : {base_pair_kind::AT, base_pair_kind::TA,
                         base_pair_kind::GC, base_pair_kind::CG})
    {
        const real_type r_min =  3.0;
        const real_type r_max = 18.0;
        const real_type dr    = (r_max - r_min) / N;

        for(std::size_t i=0; i<N; ++i)
        {
            const auto r         = r_min + i * dr;
            const auto U_attr    = potential.U_attr(bp, r);
            const auto dU_attr   = potential.dU_attr(bp, r);
            const auto U_dU_attr = potential.U_dU_attr(bp, r);

            BOOST_TEST(U_dU_attr.first  == U_attr,  boost::test_tools::tolerance(h));
            BOOST_TEST(U_dU_attr.second == dU_attr, boost::test_tools::tolerance(h));
        }
    }
    for(const auto cs : {cross_stack_kind::AA5, cross_stack_kind::AA3,
                         cross_stack_kind::AT5, cross_stack_kind::AT3,
                         cross_stack_kind::AG5, cross_stack_kind::AG3,
                         cross_stack_kind::AC5, cross_stack_kind::AC3,
                         cross_stack_kind::TA5, cross_stack_kind::TA3,
                         cross_stack_kind::TT5, cross_stack_kind::TT3,
                         cross_stack_kind::TG5, cross_stack_kind::TG3,
                         cross_stack_kind::TC5, cross_stack_kind::TC3,
                         cross_stack_kind::GA5, cross_stack_kind::GA3,
                         cross_stack_kind::GT5, cross_stack_kind::GT3,
                         cross_stack_kind::GG5, cross_stack_kind::GG3,
                         cross_stack_kind::GC5, cross_stack_kind::GC3,
                         cross_stack_kind::CA5, cross_stack_kind::CA3,
                         cross_stack_kind::CT5, cross_stack_kind::CT3,
                         cross_stack_kind::CG5, cross_stack_kind::CG3,
                         cross_stack_kind::CC5, cross_stack_kind::CC3})
    {
        const real_type r_min =  3.0;
        const real_type r_max = 18.0;
        const real_type dr    = (r_max - r_min) / N;

        for(std::size_t i=0; i<N; ++i)
        {
            const auto r         = r_min + i * dr;
            const auto U_attr    = potential.U_attr(cs, r);
            const auto dU_attr   = potential.dU_attr(cs, r);
            const auto U_dU_attr = potential.U_dU_attr(cs, r);

            BOOST_TEST(U_dU_attr.first  == U_attr,  boost::test_tools::tolerance(h));
            BOOST_TEST(U_dU_attr.second == dU_attr, boost::test_tools::tolerance(h));
        }
    }
}


BOOST_AUTO_TEST_CASE(bp_kind_3SPN2_BaseBaseInteraction)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");
    using real_type = double;
    using base_kind = mjolnir::parameter_3SPN2::base_kind;
    using base_pair_kind = mjolnir::parameter_3SPN2::base_pair_kind;

    mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type> potential({},
        typename mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>::ignore_group_type({})
        );

    BOOST_TEST(base_pair_kind::AT == potential.bp_kind(base_kind::A, base_kind::T));
    BOOST_TEST(base_pair_kind::TA == potential.bp_kind(base_kind::T, base_kind::A));
    BOOST_TEST(base_pair_kind::GC == potential.bp_kind(base_kind::G, base_kind::C));
    BOOST_TEST(base_pair_kind::CG == potential.bp_kind(base_kind::C, base_kind::G));
}
BOOST_AUTO_TEST_CASE(cs_kind_3SPN2_BaseBaseInteraction)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");
    using real_type = double;
    using base_kind = mjolnir::parameter_3SPN2::base_kind;
    using cross_stack_kind = mjolnir::parameter_3SPN2::cross_stack_kind;

    mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type> potential({},
        typename mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>::ignore_group_type({})
        );

    BOOST_TEST(cross_stack_kind::AA5 == potential.cs5_kind(base_kind::A, base_kind::A));
    BOOST_TEST(cross_stack_kind::AT5 == potential.cs5_kind(base_kind::A, base_kind::T));
    BOOST_TEST(cross_stack_kind::AG5 == potential.cs5_kind(base_kind::A, base_kind::G));
    BOOST_TEST(cross_stack_kind::AC5 == potential.cs5_kind(base_kind::A, base_kind::C));
    BOOST_TEST(cross_stack_kind::TA5 == potential.cs5_kind(base_kind::T, base_kind::A));
    BOOST_TEST(cross_stack_kind::TT5 == potential.cs5_kind(base_kind::T, base_kind::T));
    BOOST_TEST(cross_stack_kind::TG5 == potential.cs5_kind(base_kind::T, base_kind::G));
    BOOST_TEST(cross_stack_kind::TC5 == potential.cs5_kind(base_kind::T, base_kind::C));
    BOOST_TEST(cross_stack_kind::GA5 == potential.cs5_kind(base_kind::G, base_kind::A));
    BOOST_TEST(cross_stack_kind::GT5 == potential.cs5_kind(base_kind::G, base_kind::T));
    BOOST_TEST(cross_stack_kind::GG5 == potential.cs5_kind(base_kind::G, base_kind::G));
    BOOST_TEST(cross_stack_kind::GC5 == potential.cs5_kind(base_kind::G, base_kind::C));
    BOOST_TEST(cross_stack_kind::CA5 == potential.cs5_kind(base_kind::C, base_kind::A));
    BOOST_TEST(cross_stack_kind::CT5 == potential.cs5_kind(base_kind::C, base_kind::T));
    BOOST_TEST(cross_stack_kind::CG5 == potential.cs5_kind(base_kind::C, base_kind::G));
    BOOST_TEST(cross_stack_kind::CC5 == potential.cs5_kind(base_kind::C, base_kind::C));

    BOOST_TEST(cross_stack_kind::AA3 == potential.cs3_kind(base_kind::A, base_kind::A));
    BOOST_TEST(cross_stack_kind::AT3 == potential.cs3_kind(base_kind::A, base_kind::T));
    BOOST_TEST(cross_stack_kind::AG3 == potential.cs3_kind(base_kind::A, base_kind::G));
    BOOST_TEST(cross_stack_kind::AC3 == potential.cs3_kind(base_kind::A, base_kind::C));
    BOOST_TEST(cross_stack_kind::TA3 == potential.cs3_kind(base_kind::T, base_kind::A));
    BOOST_TEST(cross_stack_kind::TT3 == potential.cs3_kind(base_kind::T, base_kind::T));
    BOOST_TEST(cross_stack_kind::TG3 == potential.cs3_kind(base_kind::T, base_kind::G));
    BOOST_TEST(cross_stack_kind::TC3 == potential.cs3_kind(base_kind::T, base_kind::C));
    BOOST_TEST(cross_stack_kind::GA3 == potential.cs3_kind(base_kind::G, base_kind::A));
    BOOST_TEST(cross_stack_kind::GT3 == potential.cs3_kind(base_kind::G, base_kind::T));
    BOOST_TEST(cross_stack_kind::GG3 == potential.cs3_kind(base_kind::G, base_kind::G));
    BOOST_TEST(cross_stack_kind::GC3 == potential.cs3_kind(base_kind::G, base_kind::C));
    BOOST_TEST(cross_stack_kind::CA3 == potential.cs3_kind(base_kind::C, base_kind::A));
    BOOST_TEST(cross_stack_kind::CT3 == potential.cs3_kind(base_kind::C, base_kind::T));
    BOOST_TEST(cross_stack_kind::CG3 == potential.cs3_kind(base_kind::C, base_kind::G));
    BOOST_TEST(cross_stack_kind::CC3 == potential.cs3_kind(base_kind::C, base_kind::C));

}

