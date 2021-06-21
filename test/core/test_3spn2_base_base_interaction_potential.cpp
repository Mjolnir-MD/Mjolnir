#define BOOST_TEST_MODULE "test_3SPN2BaseBaseInteraction_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteractionPotential.hpp>

using parameters_to_test_double = std::tuple<
    mjolnir::ThreeSPN2BaseBaseGlobalPotentialParameter<double>,
    mjolnir::ThreeSPN2CBaseBaseGlobalPotentialParameter<double>
>;

BOOST_AUTO_TEST_CASE_TEMPLATE(f_3SPN2_BaseBaseInteraction_double,
        ParameterSet, parameters_to_test_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");

    using traits_type    = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type      = traits_type::real_type;
    using potential_type = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;
    using parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<traits_type>;
    using parameter_type = typename parameter_list_type::parameter_type;
    using base_kind      = mjolnir::parameter_3SPN2::base_kind  ;

    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal() / 1000.0) *
                          unit_type::avogadro_constant());
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal()) /
                            unit_type::avogadro_constant());
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom());

        phys_type::set_m_to_length(unit_type::m_to_angstrom());
        phys_type::set_length_to_m(unit_type::angstrom_to_m());

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom(), 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m(), 3));

        phys_type::set_length_unit(length);
    }

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type  pi = mjolnir::math::constants<real_type>::pi();

    constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();
    parameter_list_type parameter_list(ParameterSet{}, {
            {1, parameter_type{base_kind::A, 0, 0, invalid, invalid}},
            {3, parameter_type{base_kind::T, 1, 2, invalid, invalid}}
        }, {},
        parameter_list_type::ignore_molecule_type("Nothing"),
        parameter_list_type::ignore_group_type({})
        );

    const real_type theta_0   = 0.5  * pi;
    const real_type theta_min = 0.01 * pi;
    const real_type theta_max = 0.99 * pi;
    const real_type dtheta    = (theta_max - theta_min) / N;

    const auto K         = parameter_list.K_BP();
    const auto pi_over_K = parameter_list.pi_over_K_BP();
    const potential_type potential;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type theta = theta_min + i * dtheta;
        const real_type f1  = potential.f(K, pi_over_K, theta + h, theta_0);
        const real_type f2  = potential.f(K, pi_over_K, theta - h, theta_0);
        const real_type df_numeric  = (f1 - f2) / (2 * h);
        const real_type df_analytic = potential.df(K, pi_over_K, theta, theta_0);

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

BOOST_AUTO_TEST_CASE_TEMPLATE(U_rep_3SPN2_BaseBaseInteraction_double,
        ParameterSet, parameters_to_test_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");

    using traits_type    = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type      = traits_type::real_type;
    using potential_type = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;
    using parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<traits_type>;
    using parameter_type = typename parameter_list_type::parameter_type;
    using base_kind      = mjolnir::parameter_3SPN2::base_kind;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind;
    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal() / 1000.0) *
                          unit_type::avogadro_constant());
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal()) /
                            unit_type::avogadro_constant());
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom());

        phys_type::set_m_to_length(unit_type::m_to_angstrom());
        phys_type::set_length_to_m(unit_type::angstrom_to_m());

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom(), 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m(), 3));

        phys_type::set_length_unit(length);
    }

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();

    parameter_list_type parameter_list(ParameterSet{}, {
            {1, parameter_type{base_kind::A, 0, 0, invalid, invalid}},
            {3, parameter_type{base_kind::T, 1, 2, invalid, invalid}}
        }, {},
        parameter_list_type::ignore_molecule_type("Nothing"),
        parameter_list_type::ignore_group_type({})
        );

    const real_type r_min =  3.0;
    const real_type r_max = 12.0;
    const real_type dr    = (r_max - r_min) / N;

    const auto alpha   = parameter_list.alpha  (base_pair_kind::AT);
    const auto epsilon = parameter_list.epsilon(base_pair_kind::AT);
    const auto r0      = parameter_list.r0     (base_pair_kind::AT);

    potential_type potential;

    for(std::size_t i=0; i<N; ++i)
    {
        const auto r         = r_min + i * dr;
        const real_type pot1 = potential.U_rep(epsilon, alpha, r + h, r0);
        const real_type pot2 = potential.U_rep(epsilon, alpha, r - h, r0);
        const real_type df_numeric  = (pot1 - pot2) / (2 * h);
        const real_type df_analytic = potential.dU_rep(epsilon, alpha, r, r0);

        if(pot1 != 0.0 && pot2 == 0.0)
        {
            // the numeric differentiation becomes unstable.
        }
        else if(df_analytic == 0.0)
        {
            // derivative of repulsive part of morse potential takes zero at
            // r > r_cutoff.
            BOOST_TEST(potential.U_rep(epsilon, alpha, r, r0) == 0.0,
                       boost::test_tools::tolerance(h));
        }
        else
        {
            BOOST_TEST(df_numeric == df_analytic, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(U_attr_3SPN2_BaseBaseInteraction_double,
        ParameterSet, parameters_to_test_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");

    using traits_type    = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type      = traits_type::real_type;
    using potential_type = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;
    using parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<traits_type>;
    using parameter_type = typename parameter_list_type::parameter_type;
    using base_kind      = mjolnir::parameter_3SPN2::base_kind;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind;
    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal() / 1000.0) *
                          unit_type::avogadro_constant());
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal()) /
                            unit_type::avogadro_constant());
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom());

        phys_type::set_m_to_length(unit_type::m_to_angstrom());
        phys_type::set_length_to_m(unit_type::angstrom_to_m());

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom(), 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m(), 3));

        phys_type::set_length_unit(length);
    }

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();

    parameter_list_type parameter_list(ParameterSet{}, {
            {1, parameter_type{base_kind::A, 0, 0, invalid, invalid}},
            {3, parameter_type{base_kind::T, 1, 2, invalid, invalid}}
        }, {},
        parameter_list_type::ignore_molecule_type("Nothing"),
        parameter_list_type::ignore_group_type({})
        );

    const real_type r_min =  3.0; // avoiding the long tail around cutoff;
    const real_type r_max = 12.0; // not take numerical error around zero into account
    const real_type dr    = (r_max - r_min) / N;

    const auto alpha   = parameter_list.alpha  (base_pair_kind::AT);
    const auto epsilon = parameter_list.epsilon(base_pair_kind::AT);
    const auto r0      = parameter_list.r0     (base_pair_kind::AT);

    potential_type potential;

    for(std::size_t i=0; i<N; ++i)
    {
        const auto r         = r_min + i * dr;
        const real_type pot1 = potential.U_attr(epsilon, alpha, r + h, r0);
        const real_type pot2 = potential.U_attr(epsilon, alpha, r - h, r0);
        const real_type df_numeric  = (pot1 - pot2) / (2 * h);
        const real_type df_analytic = potential.dU_attr(epsilon, alpha, r, r0);

        if(pot1 != 0.0 && pot2 == 0.0)
        {
            // the numeric differentiation becomes unstable.
        }
        else if(df_analytic == 0.0)
        {
            // dU == 0 means that r is in `U == -epsilon` part.
            BOOST_TEST(potential.U_attr(epsilon, alpha, r, r0) == -epsilon,
                       boost::test_tools::tolerance(h));
        }
        else
        {
            BOOST_TEST(df_numeric == df_analytic, boost::test_tools::tolerance(1e-5));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(U_dU_attr_3SPN2_BaseBaseInteraction_double,
        ParameterSet, parameters_to_test_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");

    using traits_type    = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type      = traits_type::real_type;
    using potential_type = mjolnir::ThreeSPN2BaseBaseInteractionPotential<real_type>;
    using parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<traits_type>;
    using parameter_type = typename parameter_list_type::parameter_type;
    using base_kind      = mjolnir::parameter_3SPN2::base_kind;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind  ;
    using cross_stack_kind = mjolnir::parameter_3SPN2::cross_stack_kind;
    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal() / 1000.0) *
                          unit_type::avogadro_constant());
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal()) /
                            unit_type::avogadro_constant());
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom());

        phys_type::set_m_to_length(unit_type::m_to_angstrom());
        phys_type::set_length_to_m(unit_type::angstrom_to_m());

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom(), 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m(), 3));

        phys_type::set_length_unit(length);
    }

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-8;
    constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();

    parameter_list_type parameter_list(ParameterSet{}, {
            {1, parameter_type{base_kind::A, 0, 0, invalid, invalid}},
            {3, parameter_type{base_kind::T, 1, 2, invalid, invalid}}
        }, {},
        parameter_list_type::ignore_molecule_type("Nothing"),
        parameter_list_type::ignore_group_type({})
        );

    potential_type potential;

    for(const auto bp : {base_pair_kind::AT, base_pair_kind::TA,
                         base_pair_kind::GC, base_pair_kind::CG})
    {
        const real_type r_min =  3.0;
        const real_type r_max = 18.0;
        const real_type dr    = (r_max - r_min) / N;

        const auto alpha   = parameter_list.alpha  (bp);
        const auto epsilon = parameter_list.epsilon(bp);
        const auto r0      = parameter_list.r0     (bp);

        for(std::size_t i=0; i<N; ++i)
        {
            const auto r         = r_min + i * dr;
            const auto U_attr    = potential.U_attr   (epsilon, alpha, r, r0);
            const auto dU_attr   = potential.dU_attr  (epsilon, alpha, r, r0);
            const auto U_dU_attr = potential.U_dU_attr(epsilon, alpha, r, r0);

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

        const auto alpha   = parameter_list.alpha  (cs);
        const auto epsilon = parameter_list.epsilon(cs);
        const auto r0      = parameter_list.r0     (cs);

        for(std::size_t i=0; i<N; ++i)
        {
            const auto r         = r_min + i * dr;
            const auto U_attr    = potential.U_attr   (epsilon, alpha, r, r0);
            const auto dU_attr   = potential.dU_attr  (epsilon, alpha, r, r0);
            const auto U_dU_attr = potential.U_dU_attr(epsilon, alpha, r, r0);

            BOOST_TEST(U_dU_attr.first  == U_attr,  boost::test_tools::tolerance(h));
            BOOST_TEST(U_dU_attr.second == dU_attr, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(bp_kind_3SPN2_BaseBaseInteraction,
        ParameterSet, parameters_to_test_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<traits_type>;
    using base_kind        = mjolnir::parameter_3SPN2::base_kind;
    using base_pair_kind   = mjolnir::parameter_3SPN2::base_pair_kind;

    parameter_list_type parameter_list(ParameterSet{}, {}, {},
        parameter_list_type::ignore_molecule_type("Nothing"),
        parameter_list_type::ignore_group_type({})
        );

    BOOST_TEST(base_pair_kind::AT == parameter_list.bp_kind(base_kind::A, base_kind::T));
    BOOST_TEST(base_pair_kind::TA == parameter_list.bp_kind(base_kind::T, base_kind::A));
    BOOST_TEST(base_pair_kind::GC == parameter_list.bp_kind(base_kind::G, base_kind::C));
    BOOST_TEST(base_pair_kind::CG == parameter_list.bp_kind(base_kind::C, base_kind::G));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(cs_kind_3SPN2_BaseBaseInteraction,
        ParameterSet, parameters_to_test_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_3spn2_base_base_interaction_potential.log");
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using parameter_list_type = mjolnir::ThreeSPN2BaseBaseInteractionParameterList<traits_type>;
    using base_kind = mjolnir::parameter_3SPN2::base_kind;
    using cross_stack_kind = mjolnir::parameter_3SPN2::cross_stack_kind;

    parameter_list_type parameter_list(ParameterSet{}, {}, {},
        parameter_list_type::ignore_molecule_type("Nothing"),
        parameter_list_type::ignore_group_type({})
        );

    BOOST_TEST(cross_stack_kind::AA5 == parameter_list.cs5_kind(base_kind::A, base_kind::A));
    BOOST_TEST(cross_stack_kind::AT5 == parameter_list.cs5_kind(base_kind::A, base_kind::T));
    BOOST_TEST(cross_stack_kind::AG5 == parameter_list.cs5_kind(base_kind::A, base_kind::G));
    BOOST_TEST(cross_stack_kind::AC5 == parameter_list.cs5_kind(base_kind::A, base_kind::C));
    BOOST_TEST(cross_stack_kind::TA5 == parameter_list.cs5_kind(base_kind::T, base_kind::A));
    BOOST_TEST(cross_stack_kind::TT5 == parameter_list.cs5_kind(base_kind::T, base_kind::T));
    BOOST_TEST(cross_stack_kind::TG5 == parameter_list.cs5_kind(base_kind::T, base_kind::G));
    BOOST_TEST(cross_stack_kind::TC5 == parameter_list.cs5_kind(base_kind::T, base_kind::C));
    BOOST_TEST(cross_stack_kind::GA5 == parameter_list.cs5_kind(base_kind::G, base_kind::A));
    BOOST_TEST(cross_stack_kind::GT5 == parameter_list.cs5_kind(base_kind::G, base_kind::T));
    BOOST_TEST(cross_stack_kind::GG5 == parameter_list.cs5_kind(base_kind::G, base_kind::G));
    BOOST_TEST(cross_stack_kind::GC5 == parameter_list.cs5_kind(base_kind::G, base_kind::C));
    BOOST_TEST(cross_stack_kind::CA5 == parameter_list.cs5_kind(base_kind::C, base_kind::A));
    BOOST_TEST(cross_stack_kind::CT5 == parameter_list.cs5_kind(base_kind::C, base_kind::T));
    BOOST_TEST(cross_stack_kind::CG5 == parameter_list.cs5_kind(base_kind::C, base_kind::G));
    BOOST_TEST(cross_stack_kind::CC5 == parameter_list.cs5_kind(base_kind::C, base_kind::C));

    BOOST_TEST(cross_stack_kind::AA3 == parameter_list.cs3_kind(base_kind::A, base_kind::A));
    BOOST_TEST(cross_stack_kind::AT3 == parameter_list.cs3_kind(base_kind::A, base_kind::T));
    BOOST_TEST(cross_stack_kind::AG3 == parameter_list.cs3_kind(base_kind::A, base_kind::G));
    BOOST_TEST(cross_stack_kind::AC3 == parameter_list.cs3_kind(base_kind::A, base_kind::C));
    BOOST_TEST(cross_stack_kind::TA3 == parameter_list.cs3_kind(base_kind::T, base_kind::A));
    BOOST_TEST(cross_stack_kind::TT3 == parameter_list.cs3_kind(base_kind::T, base_kind::T));
    BOOST_TEST(cross_stack_kind::TG3 == parameter_list.cs3_kind(base_kind::T, base_kind::G));
    BOOST_TEST(cross_stack_kind::TC3 == parameter_list.cs3_kind(base_kind::T, base_kind::C));
    BOOST_TEST(cross_stack_kind::GA3 == parameter_list.cs3_kind(base_kind::G, base_kind::A));
    BOOST_TEST(cross_stack_kind::GT3 == parameter_list.cs3_kind(base_kind::G, base_kind::T));
    BOOST_TEST(cross_stack_kind::GG3 == parameter_list.cs3_kind(base_kind::G, base_kind::G));
    BOOST_TEST(cross_stack_kind::GC3 == parameter_list.cs3_kind(base_kind::G, base_kind::C));
    BOOST_TEST(cross_stack_kind::CA3 == parameter_list.cs3_kind(base_kind::C, base_kind::A));
    BOOST_TEST(cross_stack_kind::CT3 == parameter_list.cs3_kind(base_kind::C, base_kind::T));
    BOOST_TEST(cross_stack_kind::CG3 == parameter_list.cs3_kind(base_kind::C, base_kind::G));
    BOOST_TEST(cross_stack_kind::CC3 == parameter_list.cs3_kind(base_kind::C, base_kind::C));

}

