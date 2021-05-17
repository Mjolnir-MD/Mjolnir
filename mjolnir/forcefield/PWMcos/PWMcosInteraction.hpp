#ifndef MJOLNIR_FORCEFIELD_PWMCOS_PWMCOS_INTERACTION_HPP
#define MJOLNIR_FORCEFIELD_PWMCOS_PWMCOS_INTERACTION_HPP
#include <mjolnir/forcefield/PWMcos/PWMcosPotential.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/range.hpp>

namespace mjolnir
{

// Protein-DNA PWM-based interaction that represents sequence specificity of a
// protein.
// This is an implementation of the potential developed in the following paper.
// - C.Tan, and S.Takada (2018) JCTC
//
// U(r, theta, phi) = e(b) f(r) g(theta_1) g(theta_2) g(theta_3)
//
// e(b) represents sequence-specific interaction strength.
// f(r)     = exp(-(r-r0)^2 / 2sigma^2)
// g(theta) = 1                                  ...       |theta-theta0| <  phi
//            1 - cos^2(pi(theta-theta0) / 2phi) ... phi < |theta-theta0| < 2phi
//            0                                  ... otherwise

template<typename traitsT>
class PWMcosInteraction final : public GlobalInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using base_type       = GlobalInteractionBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using topology_type   = typename base_type::topology_type;
    using boundary_type   = typename base_type::boundary_type;
    using potential_type  = PWMcosPotential<traits_type>;
    using partition_type  = SpatialPartition<traitsT, potential_type>;

  public:
    PWMcosInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}
    ~PWMcosInteraction() {}

    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is PWMcos");

        this->potential_.initialize(sys, topol);
        this->partition_.initialize(sys, this->potential_);
        return;
    }

    void update(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is PWMcos");

        this->potential_.update(sys, topol);
        this->partition_.initialize(sys, this->potential_);
        return;
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->partition_.reduce_margin(dmargin, sys, this->potential_);
        return ;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->partition_.scale_margin(scale, sys, this->potential_);
        return ;
    }

    void      calc_force (system_type& sys)           const noexcept override;
    real_type calc_energy(const system_type& sys)     const noexcept override;
    real_type calc_force_and_energy(system_type& sys) const noexcept override;

    std::string name() const override {return "PWMcosInteraction";}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type&       potential()       noexcept {return potential_;}

    base_type* clone() const override
    {
        return new PWMcosInteraction(
                potential_type(potential_), partition_type(partition_));
    }

  private:

    potential_type potential_;
    partition_type partition_;
};

template<typename traitsT>
void PWMcosInteraction<traitsT>::calc_force(system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();
    //   DNA        protein    |
    //        i-1              | r is the length of vector `v0`.
    //     S--B       C        | angle1 is an angle formed by v0 and v1.
    //    /   |        \       | angle2 is an angle formed by v0 and v2.
    //   P    | v2      C j-1  | angle3 is an angle formed by v0 and v3.
    //    \   |        /A      |
    //  v1 S<-+-Bi--> Cj| v3   | v0: Bi -> Cj
    //    /   |    v0  \|      | v1: Bi -> Si
    //   P    |         C j+1  | v2: Bi-1 -> Bi+1
    //    \   v v2     /       | v3: Cj+1 -> Cj-1
    //     S--B       C        |
    //        i+1              |
    //
    // d/dq e(b) f(r) g(theta1) g(theta2) g(theta3)
    // = e(b) df/dr dr/dq g(theta1) g(theta2) g(theta3) +
    //   e(b) f(r) dg/dtheta1 dtheta1/dq g(theta2) g(theta3) +
    //   e(b) f(r) g(theta1) dg/dtheta2 dtheta2/dq g(theta3) +
    //   e(b) f(r) g(theta1) g(theta2) dg/dtheta3 dtheta3/dq

    constexpr auto abs_tol = math::abs_tolerance<real_type>();

    const auto energy_unit  = potential_.energy_unit();  // overall coefficient
    const auto energy_shift = potential_.energy_shift(); // overall energy shift

    for(std::size_t i=0; i < this->potential_.contacts().size(); ++i)
    {
        const auto& para = potential_.contacts()[i];
        const auto& PWM  = para.PWM;

        const auto  Ca  = para.Ca;  // C-alpha (not a calcium!)
        const auto& rCa = sys.position(Ca);

        const auto  CaN = para.CaN; // N-term Ca (Cj-1)
        const auto  CaC = para.CaC; // C-term Ca (Cj+1)

        MJOLNIR_LOG_DEBUG("Calpha = ", Ca);
        for(const auto& ptnr : this->partition_.partners(Ca))
        {
            const auto  B  = ptnr.index;          // DNA Base
            const auto  S  = ptnr.parameter().S;  // corresponding Sugar
            const auto  B5 = ptnr.parameter().B5; // Base (Bi-1)
            const auto  B3 = ptnr.parameter().B3; // Base (Bi+1)
            MJOLNIR_LOG_DEBUG("Base = ", B);

            const auto& rB     = sys.position(B);
            const auto rBCa    = sys.adjust_direction(rB, rCa); // Bi -> Cj
            const auto lBCa_sq = math::length_sq(rBCa);

            if(para.r_cut_sq < lBCa_sq) {continue;}
            MJOLNIR_LOG_DEBUG("within the cutoff");

            // ----------------------------------------------------------------
            // calculates the distance term (f(r))

            const auto rlBCa = math::rsqrt(lBCa_sq);
            const auto  lBCa = lBCa_sq * rlBCa;
            const auto  f_df = potential_.f_df(para.r0, lBCa);

            // ----------------------------------------------------------------
            // calculates the theta1 term

            const auto& rS    = sys.position(S);
            const auto  rBS   = sys.adjust_direction(rB, rS); // Bi -> Si
            const auto rlBS   = math::rlength(rBS);
            const auto dot1   = math::dot_product(rBS, rBCa);
            const auto cos1   = dot1 * rlBS * rlBCa;
            const auto theta1 = std::acos(math::clamp<real_type>(cos1,-1,1));
            const auto g_dg_1 = potential_.g_dg(para.theta1_0, theta1);

            if(g_dg_1.first == 0 && g_dg_1.second == 0)
            {
                MJOLNIR_LOG_DEBUG("theta1 - theta1_0 > 2phi");
                continue;
            }

            // ----------------------------------------------------------------
            // calculates the theta2 term

            const auto& rB3   = sys.position(B3);
            const auto& rB5   = sys.position(B5);               // (B5) -> (B3)
            const auto  rB53  = sys.adjust_direction(rB5, rB3); // Bi-1 -> Bi+1
            const auto rlB53  = math::rlength(rB53);
            const auto dot2   = math::dot_product(rB53, rBCa);
            const auto cos2   = dot2 * rlB53 * rlBCa;
            const auto theta2 = std::acos(math::clamp<real_type>(cos2,-1,1));
            const auto g_dg_2 = potential_.g_dg(para.theta2_0, theta2);

            if(g_dg_2.first == 0 && g_dg_2.second == 0)
            {
                MJOLNIR_LOG_DEBUG("theta2 - theta2_0 > 2phi");
                continue;
            }

            // ----------------------------------------------------------------
            // calculates the theta3 term

            const auto& rCaN  = sys.position(CaN);
            const auto& rCaC  = sys.position(CaC);                //(CaC) -> (CaN)
            const auto  rCCN  = sys.adjust_direction(rCaC, rCaN); // Cj+1 -> Cj-1
            const auto rlCCN  = math::rlength(rCCN);
            const auto dot3   = math::dot_product(rCCN, rBCa);
            const auto cos3   = dot3 * rlCCN * rlBCa;
            const auto theta3 = std::acos(math::clamp<real_type>(cos3,-1,1));
            const auto g_dg_3 = potential_.g_dg(para.theta3_0, theta3);

            if(g_dg_3.first == 0 && g_dg_3.second == 0)
            {
                MJOLNIR_LOG_DEBUG("theta3 - theta3_0 > 2phi");
                continue;
            }

            // ================================================================
            // calculates the force direction

            const auto Bk    = static_cast<std::size_t>(ptnr.parameter().base);
            const auto e_pwm = PWM[Bk];

            const auto coef  = para.gamma   * energy_unit;
            const auto shift = para.epsilon + energy_shift;

            const auto factor = coef * (e_pwm + shift);

            auto F_Ca  = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_B   = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_S   = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_B5  = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_B3  = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_CaN = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_CaC = math::make_coordinate<coordinate_type>(0, 0, 0);

            // ----------------------------------------------------------------
            // calculate direction term
            // = e(b) df/dr dr/dq g(theta1) g(theta2) g(theta3)

            if(g_dg_1.first != 0 && g_dg_2.first != 0 && g_dg_3.first != 0)
            {
                const auto magnitude = factor *
                    f_df.second * g_dg_1.first * g_dg_2.first * g_dg_3.first;

                F_Ca -= (magnitude * rlBCa) * rBCa;
                F_B  += (magnitude * rlBCa) * rBCa;
            }

            // ----------------------------------------------------------------
            // calculate angle1 term
            // = e(b) f(r) dg/dtheta1 dtheta1/dq g(theta2) g(theta3)

            if(g_dg_1.second != 0 && g_dg_2.first != 0 && g_dg_3.first != 0)
            {
                const auto magnitude = factor *
                    f_df.first * g_dg_1.second * g_dg_2.first * g_dg_3.first;

                const auto sin_theta = std::sin(theta1);
                const auto coef_rsin = magnitude / std::max(sin_theta, abs_tol);

                const auto Fi = (coef_rsin * rlBCa) *
                                ((cos1 * rlBCa) * rBCa - rlBS  * rBS );
                const auto Fk = (coef_rsin * rlBS ) *
                                ((cos1 * rlBS ) * rBS  - rlBCa * rBCa);

                F_Ca -= Fi;
                F_B  += (Fi + Fk);
                F_S  -= Fk;
            }

            // ----------------------------------------------------------------
            // calculate angle2 term
            // = e(b) f(r) g(theta1) dg/dtheta2 dtheta2/dq g(theta3)

            if(g_dg_1.first != 0 && g_dg_2.second != 0 && g_dg_3.first != 0)
            {
                const auto magnitude = factor *
                    f_df.first * g_dg_1.first * g_dg_2.second * g_dg_3.first;

                const auto sin_theta = std::sin(theta2);
                const auto coef_rsin = magnitude / std::max(sin_theta, abs_tol);

                const auto Fi = (coef_rsin * rlBCa) *
                                ((cos2 * rlBCa) * rBCa - rlB53 * rB53);
                const auto Fk = (coef_rsin * rlB53) *
                                ((cos2 * rlB53) * rB53 - rlBCa * rBCa);

                F_Ca -= Fi;
                F_B  += Fi;
                F_B5 += Fk;
                F_B3 -= Fk;
            }

            // ----------------------------------------------------------------
            // calculate angle3 term
            // = e(b) f(r) g(theta1) g(theta2) dg/dtheta3 dtheta3/dq

            if(g_dg_1.first != 0 && g_dg_2.first != 0 && g_dg_3.second != 0)
            {
                const auto magnitude = factor *
                    f_df.first * g_dg_1.first * g_dg_2.first * g_dg_3.second;

                const auto sin_theta = std::sin(theta3);
                const auto coef_rsin = magnitude / std::max(sin_theta, abs_tol);

                const auto Fi = (coef_rsin * rlBCa) *
                                ((cos3 * rlBCa) * rBCa - rlCCN * rCCN);
                const auto Fk = (coef_rsin * rlCCN) *
                                ((cos3 * rlCCN) * rCCN - rlBCa * rBCa);

                F_Ca  -= Fi;
                F_B   += Fi;
                F_CaC += Fk;
                F_CaN -= Fk;
            }

            // ----------------------------------------------------------------
            // collect force on Ca and B

            sys.force(Ca)  += F_Ca;
            sys.force(B )  += F_B;
            sys.force(S)   += F_S;
            sys.force(B5)  += F_B5;
            sys.force(B3)  += F_B3;
            sys.force(CaN) += F_CaN;
            sys.force(CaC) += F_CaC;

            sys.virial() += math::tensor_product(              rCa,        F_Ca)
                         +  math::tensor_product(sys.transpose(rB,   rCa), F_B)
                         +  math::tensor_product(sys.transpose(rS,   rCa), F_S)
                         +  math::tensor_product(sys.transpose(rB5,  rCa), F_B5)
                         +  math::tensor_product(sys.transpose(rB3,  rCa), F_B3)
                         +  math::tensor_product(sys.transpose(rCaN, rCa), F_CaN)
                         +  math::tensor_product(sys.transpose(rCaC, rCa), F_CaC);
        }
    }
    return ;
}

template<typename traitsT>
typename PWMcosInteraction<traitsT>::real_type
PWMcosInteraction<traitsT>::calc_energy(const system_type& sys) const noexcept
{
    const auto energy_unit  = potential_.energy_unit();  // overall coefficient
    const auto energy_shift = potential_.energy_shift(); // overall energy shift

    real_type E = 0.0;

    for(std::size_t i=0; i < this->potential_.contacts().size(); ++i)
    {
        const auto& para = potential_.contacts()[i];
        const auto& PWM  = para.PWM;

        const auto  Ca  = para.Ca;  // C-alpha (not a calcium!)
        const auto& rCa = sys.position(Ca);

        const auto  CaN = para.CaN; // N-term Ca (Cj-1)
        const auto  CaC = para.CaC; // C-term Ca (Cj+1)

        for(const auto& ptnr : this->partition_.partners(Ca))
        {
            const auto  B  = ptnr.index;          // DNA Base
            const auto  S  = ptnr.parameter().S;  // corresponding Sugar
            const auto  B5 = ptnr.parameter().B5; // Base (Bi-1)
            const auto  B3 = ptnr.parameter().B3; // Base (Bi+1)

            const auto& rB     = sys.position(B);
            const auto rBCa    = sys.adjust_direction(rB, rCa); // Bi -> Cj
            const auto lBCa_sq = math::length_sq(rBCa);

            if(para.r_cut_sq < lBCa_sq) {continue;}

            // ----------------------------------------------------------------
            // calculates the distance term (f(r))

            const auto rlBCa = math::rsqrt(lBCa_sq);
            const auto  lBCa = lBCa_sq * rlBCa;
            const auto     f = potential_.f(para.r0, lBCa);

            // ----------------------------------------------------------------
            // calculates the theta1 term

            const auto& rS    = sys.position(S);
            const auto  rBS   = sys.adjust_direction(rB, rS); // Bi -> Si
            const auto rlBS   = math::rlength(rBS);
            const auto dot1   = math::dot_product(rBS, rBCa);
            const auto cos1   = dot1 * rlBS * rlBCa;
            const auto theta1 = std::acos(math::clamp<real_type>(cos1,-1,1));
            const auto g1     = potential_.g(para.theta1_0, theta1);

            if(g1 == 0) {continue;}

            // ----------------------------------------------------------------
            // calculates the theta2 term

            const auto& rB3   = sys.position(B3);
            const auto& rB5   = sys.position(B5);               // (B5) -> (B3)
            const auto  rB53  = sys.adjust_direction(rB5, rB3); // Bi-1 -> Bi+1
            const auto rlB53  = math::rlength(rB53);
            const auto dot2   = math::dot_product(rB53, rBCa);
            const auto cos2   = dot2 * rlB53 * rlBCa;
            const auto theta2 = std::acos(math::clamp<real_type>(cos2,-1,1));
            const auto g2     = potential_.g(para.theta2_0, theta2);

            if(g2 == 0) {continue;}

            // ----------------------------------------------------------------
            // calculates the theta3 term

            const auto& rCaN  = sys.position(CaN);
            const auto& rCaC  = sys.position(CaC);                //(CaC) -> (CaN)
            const auto  rCCN  = sys.adjust_direction(rCaC, rCaN); // Cj+1 -> Cj-1
            const auto rlCCN  = math::rlength(rCCN);
            const auto dot3   = math::dot_product(rCCN, rBCa);
            const auto cos3   = dot3 * rlCCN * rlBCa;
            const auto theta3 = std::acos(math::clamp<real_type>(cos3,-1,1));
            const auto g3     = potential_.g(para.theta3_0, theta3);

            if(g3 == 0) {continue;}

            // ================================================================
            // calculates the force direction

            const auto Bk    = static_cast<std::size_t>(ptnr.parameter().base);
            const auto e_pwm = PWM[Bk];
            const auto coef  = para.gamma   * energy_unit;
            const auto shift = para.epsilon + energy_shift;

            E += coef * (e_pwm + shift) * f * g1 * g2 * g3;
        }
    }
    return E;
}

template<typename traitsT>
typename PWMcosInteraction<traitsT>::real_type
PWMcosInteraction<traitsT>::calc_force_and_energy(system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();
    //   DNA        protein    |
    //        i-1              | r is the length of vector `v0`.
    //     S--B       C        | angle1 is an angle formed by v0 and v1.
    //    /   |        \       | angle2 is an angle formed by v0 and v2.
    //   P    | v2      C j-1  | angle3 is an angle formed by v0 and v3.
    //    \   |        /A      |
    //  v1 S<-+-Bi--> Cj| v3   | v0: Bi -> Cj
    //    /   |    v0  \|      | v1: Bi -> Si
    //   P    |         C j+1  | v2: Bi-1 -> Bi+1
    //    \   v v2     /       | v3: Cj+1 -> Cj-1
    //     S--B       C        |
    //        i+1              |
    //
    // d/dq e(b) f(r) g(theta1) g(theta2) g(theta3)
    // = e(b) df/dr dr/dq g(theta1) g(theta2) g(theta3) +
    //   e(b) f(r) dg/dtheta1 dtheta1/dq g(theta2) g(theta3) +
    //   e(b) f(r) g(theta1) dg/dtheta2 dtheta2/dq g(theta3) +
    //   e(b) f(r) g(theta1) g(theta2) dg/dtheta3 dtheta3/dq

    constexpr auto abs_tol = math::abs_tolerance<real_type>();

    const auto energy_unit  = potential_.energy_unit();  // overall coefficient
    const auto energy_shift = potential_.energy_shift(); // overall energy shift

    real_type energy = 0;
    for(std::size_t i=0; i < this->potential_.contacts().size(); ++i)
    {
        const auto& para = potential_.contacts()[i];
        const auto& PWM  = para.PWM;

        const auto  Ca  = para.Ca;  // C-alpha (not a calcium!)
        const auto& rCa = sys.position(Ca);

        const auto  CaN = para.CaN; // N-term Ca (Cj-1)
        const auto  CaC = para.CaC; // C-term Ca (Cj+1)

        MJOLNIR_LOG_DEBUG("Calpha = ", Ca);
        for(const auto& ptnr : this->partition_.partners(Ca))
        {
            const auto  B  = ptnr.index;          // DNA Base
            const auto  S  = ptnr.parameter().S;  // corresponding Sugar
            const auto  B5 = ptnr.parameter().B5; // Base (Bi-1)
            const auto  B3 = ptnr.parameter().B3; // Base (Bi+1)
            MJOLNIR_LOG_DEBUG("Base = ", B);

            const auto& rB     = sys.position(B);
            const auto rBCa    = sys.adjust_direction(rB, rCa); // Bi -> Cj
            const auto lBCa_sq = math::length_sq(rBCa);

            if(para.r_cut_sq < lBCa_sq) {continue;}
            MJOLNIR_LOG_DEBUG("within the cutoff");

            // ----------------------------------------------------------------
            // calculates the distance term (f(r))

            const auto rlBCa = math::rsqrt(lBCa_sq);
            const auto  lBCa = lBCa_sq * rlBCa;
            const auto  f_df = potential_.f_df(para.r0, lBCa);

            // ----------------------------------------------------------------
            // calculates the theta1 term

            const auto& rS    = sys.position(S);
            const auto  rBS   = sys.adjust_direction(rB, rS); // Bi -> Si
            const auto rlBS   = math::rlength(rBS);
            const auto dot1   = math::dot_product(rBS, rBCa);
            const auto cos1   = dot1 * rlBS * rlBCa;
            const auto theta1 = std::acos(math::clamp<real_type>(cos1,-1,1));
            const auto g_dg_1 = potential_.g_dg(para.theta1_0, theta1);

            if(g_dg_1.first == 0 && g_dg_1.second == 0)
            {
                MJOLNIR_LOG_DEBUG("theta1 - theta1_0 > 2phi");
                continue;
            }

            // ----------------------------------------------------------------
            // calculates the theta2 term

            const auto& rB3   = sys.position(B3);
            const auto& rB5   = sys.position(B5);               // (B5) -> (B3)
            const auto  rB53  = sys.adjust_direction(rB5, rB3); // Bi-1 -> Bi+1
            const auto rlB53  = math::rlength(rB53);
            const auto dot2   = math::dot_product(rB53, rBCa);
            const auto cos2   = dot2 * rlB53 * rlBCa;
            const auto theta2 = std::acos(math::clamp<real_type>(cos2,-1,1));
            const auto g_dg_2 = potential_.g_dg(para.theta2_0, theta2);

            if(g_dg_2.first == 0 && g_dg_2.second == 0)
            {
                MJOLNIR_LOG_DEBUG("theta2 - theta2_0 > 2phi");
                continue;
            }

            // ----------------------------------------------------------------
            // calculates the theta3 term

            const auto& rCaN  = sys.position(CaN);
            const auto& rCaC  = sys.position(CaC);                //(CaC) -> (CaN)
            const auto  rCCN  = sys.adjust_direction(rCaC, rCaN); // Cj+1 -> Cj-1
            const auto rlCCN  = math::rlength(rCCN);
            const auto dot3   = math::dot_product(rCCN, rBCa);
            const auto cos3   = dot3 * rlCCN * rlBCa;
            const auto theta3 = std::acos(math::clamp<real_type>(cos3,-1,1));
            const auto g_dg_3 = potential_.g_dg(para.theta3_0, theta3);

            if(g_dg_3.first == 0 && g_dg_3.second == 0)
            {
                MJOLNIR_LOG_DEBUG("theta3 - theta3_0 > 2phi");
                continue;
            }

            // ================================================================
            // calculates the force direction

            const auto Bk    = static_cast<std::size_t>(ptnr.parameter().base);
            const auto e_pwm = PWM[Bk];

            const auto coef  = para.gamma   * energy_unit;
            const auto shift = para.epsilon + energy_shift;

            const auto factor = coef * (e_pwm + shift);

            energy += factor * f_df.first *
                      g_dg_1.first * g_dg_2.first * g_dg_3.first;

            auto F_Ca  = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_B   = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_S   = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_B5  = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_B3  = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_CaN = math::make_coordinate<coordinate_type>(0, 0, 0);
            auto F_CaC = math::make_coordinate<coordinate_type>(0, 0, 0);

            // ----------------------------------------------------------------
            // calculate direction term
            // = e(b) df/dr dr/dq g(theta1) g(theta2) g(theta3)

            if(g_dg_1.first != 0 && g_dg_2.first != 0 && g_dg_3.first != 0)
            {
                const auto magnitude = factor *
                    f_df.second * g_dg_1.first * g_dg_2.first * g_dg_3.first;

                F_Ca -= (magnitude * rlBCa) * rBCa;
                F_B  += (magnitude * rlBCa) * rBCa;
            }

            // ----------------------------------------------------------------
            // calculate angle1 term
            // = e(b) f(r) dg/dtheta1 dtheta1/dq g(theta2) g(theta3)

            if(g_dg_1.second != 0 && g_dg_2.first != 0 && g_dg_3.first != 0)
            {
                const auto magnitude = factor *
                    f_df.first * g_dg_1.second * g_dg_2.first * g_dg_3.first;

                const auto sin_theta = std::sin(theta1);
                const auto coef_rsin = magnitude / std::max(sin_theta, abs_tol);

                const auto Fi = (coef_rsin * rlBCa) *
                                ((cos1 * rlBCa) * rBCa - rlBS  * rBS );
                const auto Fk = (coef_rsin * rlBS ) *
                                ((cos1 * rlBS ) * rBS  - rlBCa * rBCa);

                F_Ca -= Fi;
                F_B  += (Fi + Fk);
                F_S  -= Fk;
            }

            // ----------------------------------------------------------------
            // calculate angle2 term
            // = e(b) f(r) g(theta1) dg/dtheta2 dtheta2/dq g(theta3)

            if(g_dg_1.first != 0 && g_dg_2.second != 0 && g_dg_3.first != 0)
            {
                const auto magnitude = factor *
                    f_df.first * g_dg_1.first * g_dg_2.second * g_dg_3.first;

                const auto sin_theta = std::sin(theta2);
                const auto coef_rsin = magnitude / std::max(sin_theta, abs_tol);

                const auto Fi = (coef_rsin * rlBCa) *
                                ((cos2 * rlBCa) * rBCa - rlB53 * rB53);
                const auto Fk = (coef_rsin * rlB53) *
                                ((cos2 * rlB53) * rB53 - rlBCa * rBCa);

                F_Ca -= Fi;
                F_B  += Fi;
                F_B5 += Fk;
                F_B3 -= Fk;
            }

            // ----------------------------------------------------------------
            // calculate angle3 term
            // = e(b) f(r) g(theta1) g(theta2) dg/dtheta3 dtheta3/dq

            if(g_dg_1.first != 0 && g_dg_2.first != 0 && g_dg_3.second != 0)
            {
                const auto magnitude = factor *
                    f_df.first * g_dg_1.first * g_dg_2.first * g_dg_3.second;

                const auto sin_theta = std::sin(theta3);
                const auto coef_rsin = magnitude / std::max(sin_theta, abs_tol);

                const auto Fi = (coef_rsin * rlBCa) *
                                ((cos3 * rlBCa) * rBCa - rlCCN * rCCN);
                const auto Fk = (coef_rsin * rlCCN) *
                                ((cos3 * rlCCN) * rCCN - rlBCa * rBCa);

                F_Ca  -= Fi;
                F_B   += Fi;
                F_CaC += Fk;
                F_CaN -= Fk;
            }

            // ----------------------------------------------------------------
            // collect force on Ca and B

            sys.force(Ca)  += F_Ca;
            sys.force(B )  += F_B ;
            sys.force(S)   += F_S;
            sys.force(B5)  += F_B5;
            sys.force(B3)  += F_B3;
            sys.force(CaN) += F_CaN;
            sys.force(CaC) += F_CaC;

            sys.virial() += math::tensor_product(              rCa,        F_Ca)
                         +  math::tensor_product(sys.transpose(rB,   rCa), F_B)
                         +  math::tensor_product(sys.transpose(rS,   rCa), F_S)
                         +  math::tensor_product(sys.transpose(rB5,  rCa), F_B5)
                         +  math::tensor_product(sys.transpose(rB3,  rCa), F_B3)
                         +  math::tensor_product(sys.transpose(rCaN, rCa), F_CaN)
                         +  math::tensor_product(sys.transpose(rCaC, rCa), F_CaC);
        }
    }
    return energy;
}


} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{

extern template class PWMcosInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class PWMcosInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class PWMcosInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class PWMcosInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif// MJOLNIR_FORCEFIELD_PWMCOS_PWMCOS_INTERACTION_HPP
