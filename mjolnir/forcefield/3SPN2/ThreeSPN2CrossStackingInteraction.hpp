#ifndef MJOLNIR_FORCEFIELD_3SPN2_CROSS_STACKING_INTERACTION_HPP
#define MJOLNIR_FORCEFIELD_3SPN2_CROSS_STACKING_INTERACTION_HPP
#include <mjolnir/forcefield/3SPN2/ThreeSPN2CrossStackingPotential.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <memory>

namespace mjolnir
{

// 3SPN.2 Base-Base non-local interaction (base pairing & cross-stacking).
// This is an implementation of the potential developed in the following paper.
// - D.M.Hinckley, G.S.Freeman, J.K.Whitmer, and J.J.de Pablo (2013) J. Chem. Phys.
//   doi: 10.1063/1.4822042
//
// This interaction is deeply coupled with its potential, so it does not receive
// a potential type as a template argument. It always uses the same potential,
// `ThreeSPN2CrossStackingInteractionPotential`.
//
// It shares CellList between BasePairing and CrossStacking.
// It constructs CellList for BasePairing and re-use the pairs for CrossStacking.
//
template<typename traitsT>
class ThreeSPN2CrossStackingInteraction final : public GlobalInteractionBase<traitsT>
{
  public:

    using traits_type         = traitsT;
    using base_type           = GlobalInteractionBase<traits_type>;
    using real_type           = typename base_type::real_type;
    using coordinate_type     = typename base_type::coordinate_type;
    using system_type         = typename base_type::system_type;
    using topology_type       = typename base_type::topology_type;
    using boundary_type       = typename base_type::boundary_type;
    using potential_type      = ThreeSPN2CrossStackingPotential<real_type>;
    using parameter_list_type = ThreeSPN2CrossStackingParameterList<traits_type>;
    using partition_type      = SpatialPartition<traitsT, potential_type>;

    using base_kind        = parameter_3SPN2::bead_kind;
    using cross_stack_kind = parameter_3SPN2::cross_stack_kind;

  public:

    ThreeSPN2CrossStackingInteraction(potential_type&& pot, parameter_list_type&& para, partition_type&& part)
        : potential_(std::move(pot)), parameters_(std::move(para)), partition_(std::move(part))
    {}
    ~ThreeSPN2CrossStackingInteraction() {}

    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());

        this->potential_ .initialize(sys);
        this->parameters_.initialize(sys, topol, potential_);
        this->partition_ .initialize(sys, this->parameters_);
    }

    void update(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());

        this->potential_ .update(sys);
        this->parameters_.update(sys, topol, potential_);
        this->partition_ .initialize(sys, this->parameters_);
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->partition_.reduce_margin(dmargin, sys, this->parameters_);
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->partition_.scale_margin(scale, sys, this->parameters_);
        return;
    }

    void calc_force(system_type& sys)           const noexcept override
    {
        this->template calc_force_energy_virial_impl<false, false>(sys);
        return ;
    }
    void calc_force_and_virial(system_type& sys) const noexcept override
    {
        this->template calc_force_energy_virial_impl<false, true>(sys);
        return;
    }
    real_type calc_energy(const system_type& sys)     const noexcept override;
    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_energy_virial_impl<true, false>(sys);
    }
    real_type calc_force_virial_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_energy_virial_impl<true, true>(sys);
    }
    std::string name() const override {return "3SPN2CrossStacking";}

    parameter_list_type const& parameters() const noexcept {return parameters_;}
    parameter_list_type&       parameters()       noexcept {return parameters_;}

    base_type* clone() const override
    {
        return new ThreeSPN2CrossStackingInteraction(*this);
    }

  private:

    template<bool NeedEnergy, bool NeedVirial>
    real_type calc_force_energy_virial_impl(system_type& sys) const noexcept;

  private:

    potential_type      potential_;
    parameter_list_type parameters_;
    partition_type      partition_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

// ===========================================================================

template<typename traitsT>
typename ThreeSPN2CrossStackingInteraction<traitsT>::real_type
ThreeSPN2CrossStackingInteraction<traitsT>::calc_energy(
        const system_type& sys) const noexcept
{
    const auto K_CS            = parameters_.K_CS();
    const auto pi_over_K_CS    = parameters_.pi_over_K_CS();
    const auto K_BP            = parameters_.K_BP();
    const auto pi_over_K_BP    = parameters_.pi_over_K_BP();

    real_type E_CS = 0.0;

    for(const std::size_t Bi : this->parameters_.leading_participants())
    {
        // ================================================================
        // cross stacking
        // f(theta_3) f(theta_CS) U_attr(epsilon, alpha, rij)
        //
        //           Bi5   Bj3
        //  5'    o -- o===o -- o     3'
        //  ^    /      \ /      \    |
        //  | P o        X        o P |
        //  |    \   Bi / \      /    |
        //  |  Si o -- o===o -- o Sj  |
        //  |    /      \ / Bj   \    |
        //  | P o        X        o P |
        //  |    \      / \      /    v
        //  3'    o -- o===o -- o     5'
        //           Bi3    Bj5

        const auto& rBi = sys.position(Bi);

        const auto& para_Bi = parameters_.at(Bi);
        const auto Bi5 = para_Bi.B5;
        const auto Bi3 = para_Bi.B3;
        const auto Si  = para_Bi.S;
        const auto Bi_base = para_Bi.base;

        const bool Bi5_exists = (Bi5 != potential_type::invalid());
        const bool Bi3_exists = (Bi3 != potential_type::invalid());

        for(const auto& ptnr : this->partition_.partners(Bi))
        {
            const auto  Bj   = ptnr.index;
            const auto& rBj  = sys.position(Bj);
            const auto& para = ptnr.parameter();
            const auto  bp_kind = para.bp_kind;

            const auto  Bj_base = this->parameters_.at(Bj).base;

            const auto  Bij    = sys.adjust_direction(rBi, rBj); // Bi -> Bj
            const auto lBij_sq = math::length_sq(Bij);
            if(lBij_sq > parameters_.cutoff_sq())
            {
                continue; // if the base pair does not form a pair, then do nothing.
            }

            const auto Bj5 = parameters_.at(Bj).B5;
            const auto Bj3 = parameters_.at(Bj).B3;
            const auto Sj  = parameters_.at(Bj).S;
            const bool Bj5_exists = (Bj5 != potential_type::invalid());
            const bool Bj3_exists = (Bj3 != potential_type::invalid());

            if(!Bi3_exists && !Bi5_exists && !Bj3_exists && !Bj5_exists)
            {
                continue; // if no interacting pair exist, do nothing.
            }

            const auto& rSi = sys.position(Si);
            const auto& rSj = sys.position(Sj);
            const auto SBi     = sys.adjust_direction(rSi, rBi); // Si -> Bi
            const auto SBj     = sys.adjust_direction(rSj, rBj); // Sj -> Bj
            const auto lSBi_sq = math::length_sq(SBi); // |SBi|^2
            const auto lSBj_sq = math::length_sq(SBj); // |SBj|^2
            const auto rlSBi   = math::rsqrt(lSBi_sq); // 1 / |SBi|
            const auto rlSBj   = math::rsqrt(lSBj_sq); // 1 / |SBj|

            const auto dot_SBiSBj = math::dot_product(SBi, SBj);
            const auto cos_theta3 = dot_SBiSBj * rlSBi * rlSBj;
            const auto theta3     = std::acos(math::clamp<real_type>(cos_theta3, -1, 1));
            const auto theta3_0   = parameters_.theta3_0(bp_kind);
            const auto f3         = potential_.f(K_BP, pi_over_K_BP, theta3, theta3_0);

            if(f3 == real_type(0))
            {
                continue; // no cross-stacking is formed.
            }

            // ----------------------------------------------------------------

            const auto calc_energy_element = [this, &sys, f3, K_CS, pi_over_K_CS](
                    const cross_stack_kind cs_kind,
                    const coordinate_type& rB1, const coordinate_type& rB2,
                    const coordinate_type& SB1, const real_type rlSB1
                    ) noexcept -> real_type
            {
                //       S1   B1   Bj   Sj
                //  5'    o--> o===o <--o     3'
                //  ^    /   `--\        \    |
                //  | P o   tCS  \        o P |
                //  |    \        \      /    v
                //  3'    o -- o===o -- o     5'
                //                 B2

                const auto epsilon_cs = this->parameters_.epsilon(cs_kind);
                const auto alpha_cs   = this->parameters_.alpha(cs_kind);
                const auto r0_cs      = this->parameters_.r0(cs_kind);

                const auto B21     = sys.adjust_direction(rB2, rB1);
                const auto lB21_sq = math::length_sq(B21);
                const auto rlB21   = math::rsqrt(lB21_sq);

                const auto cos_theta_CS = math::dot_product(SB1, B21) * rlSB1 * rlB21;
                const auto theta_CS     = std::acos(
                        math::clamp<real_type>(cos_theta_CS, -1, 1));
                const auto theta_CS_0   = parameters_.thetaCS_0(cs_kind);

                const auto fCS = potential_.f(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                if(fCS == real_type(0))
                {
                    return real_type(0);
                }
                const auto lB21  = lB21_sq * rlB21;
                const auto U_attr = this->potential_.U_attr(epsilon_cs, alpha_cs, lB21, r0_cs);
                return 0.5 * f3 * fCS * U_attr;
            };

            // Si--Bi -> Bj5
            if(Bj5_exists)
            {
                //       Si   Bi   Bj   Sj
                //  5'    o--> o===o <--o     3'
                //  ^    /   `--\        \    |
                //  | P o   tCS  \        o P |
                //  |    \        \      /    v
                //  3'    o -- o===o -- o     5'
                //           Bi3   Bj5
                //
                // Bi(5' side) -> Bj5(5' side): cs5_kind

                const auto Bj5_base = parameters_.at(Bj5).base;
                const auto cs_kind  = parameters_.cs5_kind(Bi_base, Bj5_base);
                const auto& rBj5 = sys.position(Bj5);

                E_CS += calc_energy_element(cs_kind, rBi, rBj5, SBi, rlSBi);
            }
            // Sj--Bj -> Bi5
            if(Bi5_exists)
            {
                //           Bi5
                //  5'    o -- o===o -- o     3'
                //  ^    /      \        \    |
                //  | P o        \        o P |
                //  |    \        \--.   /    |
                //  |  Si o -- o===o -- o Sj  v
                //  3'       Bi     Bj        5'
                //
                // Bj(5' side) -> Bj5(5' side): cs5_kind

                const auto Bi5_base = parameters_.at(Bi5).base;
                const auto cs_kind  = parameters_.cs5_kind(Bj_base, Bi5_base);
                const auto& rBi5 = sys.position(Bi5);

                E_CS += calc_energy_element(cs_kind, rBj, rBi5, SBj, rlSBj);
            }

            if(Bi3_exists)
            {
                //       Si   Bi   Bj   Sj
                //  5'    o--> o===o <--o     3'
                //  ^    /        /_.'   \    |
                //  | P o        /  tCS   o P |
                //  |    \      /        /    v
                //  3'    o -- o===o -- o     5'
                //           Bi3   Bj5

                const auto Bi3_base = parameters_.at(Bi3).base;
                const auto cs_kind  = parameters_.cs3_kind(Bj_base, Bi3_base);
                const auto& rBi3 = sys.position(Bi3);

                E_CS += calc_energy_element(cs_kind, rBj, rBi3, SBj, rlSBj);
            }
            if(Bj3_exists)
            {
                // ----------------------------------------------------------------
                //
                //           Bi5   Bj3
                //  5'    o--> o===o <--o     3'
                //  ^    /        /_.'   \    |
                //  | P o        /  tCS   o P |
                //  |    \      /        /    v
                //  3'    o -- o===o -- o     5'
                //       Si   Bi   Bj   Sj

                const auto Bj3_base = parameters_.at(Bj3).base;
                const auto cs_kind  = parameters_.cs3_kind(Bi_base, Bj3_base);
                const auto& rBj3 = sys.position(Bj3);

                E_CS += calc_energy_element(cs_kind, rBi, rBj3, SBi, rlSBi);
            }
        }
    }
    return E_CS;
}

template<typename traitsT>
template<bool NeedEnergy, bool NeedVirial>
typename ThreeSPN2CrossStackingInteraction<traitsT>::real_type
ThreeSPN2CrossStackingInteraction<traitsT>::calc_force_energy_virial_impl(system_type& sys) const noexcept
{
    constexpr auto tolerance = math::abs_tolerance<real_type>();

    const auto K_CS            = parameters_.K_CS();
    const auto pi_over_K_CS    = parameters_.pi_over_K_CS();
    const auto K_BP            = parameters_.K_BP();
    const auto pi_over_K_BP    = parameters_.pi_over_K_BP();

    real_type energy = 0;
    for(const std::size_t Bi : this->parameters_.leading_participants())
    {
        // ================================================================
        // cross stacking
        // f(theta_3) f(theta_CS) U_attr(epsilon, alpha, rij)
        //
        //           Bi5   Bj3
        //  5'    o -- o===o -- o     3'
        //  ^    /      \ /      \    |
        //  | P o        X        o P |
        //  |    \   Bi / \      /    |
        //  |  Si o -- o===o -- o Sj  |
        //  |    /      \ / Bj   \    |
        //  | P o        X        o P |
        //  |    \      / \      /    v
        //  3'    o -- o===o -- o     5'
        //           Bi3    Bj5
        //
        // d/dr Ucs =
        //    df/dtheta3 f(theta_CS)  U_attr(eps, alp, rij) dtheta_3  /dr
        //  + f(theta_3) df/dtheta_CS U_attr(eps, alp, rij) dtheta_CS /dr
        //  + f(theta_3) f(theta_CS)  dU_attr/drij          drij/dr
        //

        const auto& rBi = sys.position(Bi);

        const auto& para_Bi = parameters_.at(Bi);
        const auto Bi5 = para_Bi.B5;
        const auto Bi3 = para_Bi.B3;
        const auto Si  = para_Bi.S;
        const auto Bi_base = para_Bi.base;

        const bool Bi5_exists = (Bi5 != potential_type::invalid());
        const bool Bi3_exists = (Bi3 != potential_type::invalid());

        for(const auto& ptnr : this->partition_.partners(Bi))
        {
            const auto  Bj   = ptnr.index;
            const auto& rBj  = sys.position(Bj);
            const auto& para = ptnr.parameter();
            const auto  bp_kind = para.bp_kind;

            const auto Bij = sys.adjust_direction(rBi, rBj); // Bi -> Bj
            const auto lBij_sq = math::length_sq(Bij);
            if(lBij_sq > parameters_.cutoff_sq())
            {
                continue;
            }

            const auto Bj5 = parameters_.at(Bj).B5;
            const auto Bj3 = parameters_.at(Bj).B3;
            const auto Sj  = parameters_.at(Bj).S;
            const bool Bj5_exists = (Bj5 != potential_type::invalid());
            const bool Bj3_exists = (Bj3 != potential_type::invalid());

            if(!Bi3_exists && !Bi5_exists && !Bj3_exists && !Bj5_exists)
            {
                continue; // if no interacting pair exist, do nothing.
            }

            const auto  Bj_base = this->parameters_.at(Bj).base;

            const auto& rSi = sys.position(Si);
            const auto& rSj = sys.position(Sj);
            const auto SBi = sys.adjust_direction(rSi, rBi); // Si -> Bi
            const auto SBj = sys.adjust_direction(rSj, rBj); // Sj -> Bj
            const auto lSBi_sq = math::length_sq(SBi); // |SBi|^2
            const auto lSBj_sq = math::length_sq(SBj); // |SBj|^2
            const auto rlSBi   = math::rsqrt(lSBi_sq); // 1 / |SBi|
            const auto rlSBj   = math::rsqrt(lSBj_sq); // 1 / |SBj|
            const auto BSi_reg = -rlSBi * SBi;
            const auto BSj_reg = -rlSBj * SBj;

            // ----------------------------------------------------------------
            // calc common part, theta3 and dtheta3/dr.

            const auto cos_theta3 = math::dot_product(BSi_reg, BSj_reg);
            const auto theta3     = std::acos(math::clamp<real_type>(cos_theta3, -1, 1));
            const auto theta3_0   = parameters_.theta3_0(bp_kind);
            const auto f3         = potential_.f(K_BP, pi_over_K_BP, theta3, theta3_0);
            if(f3 == real_type(0))
            {
                // f(theta) == 0 means df(theta) is also zero.
                // so here, both cross-stacking becomes zero. skip them.
                continue;
            }
            const auto df3 = potential_.df(K_BP, pi_over_K_BP, theta3, theta3_0);

            // ----------------------------------------------------------------
            // force directions for dtheta3/dr

            // here, theta3 is a return value of acos, so it is in [0, pi].
            // and inside it, sin(theta3) should be larger than 0.
            const auto sin_theta3  = std::sin(theta3);
            const auto rsin_theta3 = sin_theta3 > tolerance ?
                real_type(1) / sin_theta3 : real_type(1) / tolerance;

            const auto fBi_theta3 = rsin_theta3 * rlSBi * (BSj_reg - cos_theta3 * BSi_reg);
            const auto fBj_theta3 = rsin_theta3 * rlSBj * (BSi_reg - cos_theta3 * BSj_reg);
            const auto fSi_theta3 = real_type(-1) * fBi_theta3;
            const auto fSj_theta3 = real_type(-1) * fBj_theta3;

            auto f_Si  = math::make_coordinate<coordinate_type>(0,0,0);
            auto f_Sj  = math::make_coordinate<coordinate_type>(0,0,0);
            auto f_Bi  = math::make_coordinate<coordinate_type>(0,0,0);
            auto f_Bi3 = math::make_coordinate<coordinate_type>(0,0,0);
            auto f_Bi5 = math::make_coordinate<coordinate_type>(0,0,0);
            auto f_Bj  = math::make_coordinate<coordinate_type>(0,0,0);
            auto f_Bj3 = math::make_coordinate<coordinate_type>(0,0,0);
            auto f_Bj5 = math::make_coordinate<coordinate_type>(0,0,0);

            const auto calc_force_element = [this, &sys, f3, df3, K_CS, pi_over_K_CS] (
                    const cross_stack_kind cs_kind,
                    const coordinate_type& rB1, const coordinate_type& rB_next,
                    const coordinate_type& SB1, const real_type rlSB1,
                    const coordinate_type& BS1_reg,
                    // force on S1, B1, B2, S2 from theta3 term
                    const coordinate_type& fS1_theta3,
                    const coordinate_type& fB1_theta3,
                    const coordinate_type& fS2_theta3,
                    const coordinate_type& fB2_theta3,
                    // force output
                    coordinate_type& f_S1, coordinate_type& f_B1,
                    coordinate_type& f_S2, coordinate_type& f_B2,
                    coordinate_type& f_Bnext) -> real_type // returns energy
            {
                //       S1   B1   B2   S2
                //  5'    o--> o===o <--o     3'
                //  ^    /   `--\        \    |
                //  | P o   tCS  \        o P |
                //  |    \        \      /    v
                //  3'    o -- o===o -- o     5'
                //                 B_next

                const auto epsilon_cs = parameters_.epsilon(cs_kind);
                const auto alpha_cs   = parameters_.alpha(cs_kind);
                const auto r0_cs      = parameters_.r0(cs_kind);

                const auto Bn1     = sys.adjust_direction(rB_next, rB1);
                const auto lBn1_sq = math::length_sq(Bn1);
                const auto rlBn1   = math::rsqrt(lBn1_sq);

                const auto dot_theta_CS = math::dot_product(SB1, Bn1);
                const auto cos_theta_CS = dot_theta_CS * rlSB1 * rlBn1;
                const auto theta_CS     = std::acos(math::clamp<real_type>(cos_theta_CS, -1, 1));
                const auto theta_CS_0   = parameters_.thetaCS_0(cs_kind);

                const auto fCS  = potential_.f(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                const auto dfCS = potential_.df(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                if(fCS == real_type(0))
                {
                    return real_type(0);
                }

                const auto lBn1 = lBn1_sq * rlBn1;
                const auto U_dU_attr = potential_.U_dU_attr(epsilon_cs, alpha_cs, lBn1, r0_cs);
                const auto Bn1_reg  = rlBn1 * Bn1;

                // --------------------------------------------------------
                // df/dtheta3 f(theta_CS)  U_attr(eps, alp, rij) dtheta_3 /dr
                if(df3 != real_type(0))
                {
                    const auto coef = -0.5 * df3 * fCS * U_dU_attr.first;
                    f_S1 += coef * fS1_theta3;
                    f_S2 += coef * fS2_theta3;
                    f_B1 += coef * fB1_theta3;
                    f_B2 += coef * fB2_theta3;
                }
                // --------------------------------------------------------
                // f(theta_3) df/dtheta_CS U_attr(eps, alp, rij) dtheta_CS/dr
                if(dfCS != real_type(0))
                {
                    const auto coef         = -0.5 * f3 * dfCS * U_dU_attr.first;
                    const auto sin_theta_CS = std::sin(theta_CS);
                    const auto coef_rsin    = (sin_theta_CS > tolerance) ?
                               (coef / sin_theta_CS) : (coef / tolerance);

                    const auto fS1     =  coef_rsin * rlSB1 * (cos_theta_CS * BS1_reg + Bn1_reg);
                    const auto fB_next = -coef_rsin * rlBn1 * (cos_theta_CS * Bn1_reg + BS1_reg);

                    f_S1    +=  fS1;
                    f_B1    -= (fS1 + fB_next);
                    f_Bnext +=  fB_next;
                }
                // --------------------------------------------------------
                // f(theta_3) f(theta_CS)  dU_attr/drij          drij/dr
                if(U_dU_attr.second != real_type(0.0))
                {
                    const auto coef = -0.5 * f3 * fCS * U_dU_attr.second;
                    f_B1    += coef * Bn1_reg;
                    f_Bnext -= coef * Bn1_reg;
                }

                // -------------------------------------------------------------

                if /*constexpr*/ (NeedEnergy)
                {
                    return 0.5 * f3 * fCS * U_dU_attr.first;
                }
                else
                {
                    return real_type(0);
                }
            };

            // -----------------------------------------------------------------

            // Si--Bi -> Bj5
            if(Bj5_exists)
            {
                //       Si   Bi   Bj   Sj
                //  5'    o--> o===o <--o     3'
                //  ^    /   `--\        \    |
                //  | P o   tCS  \        o P |
                //  |    \        \      /    v
                //  3'    o -- o===o -- o     5'
                //           Bi3   Bj5
                //
                // Bi(5' side) -> Bj5(5' side): cs5_kind

                const auto Bj5_base = parameters_.at(Bj5).base;
                const auto cs_kind  = parameters_.cs5_kind(Bi_base, Bj5_base);
                const auto& rBj5 = sys.position(Bj5);

                energy += calc_force_element(cs_kind, rBi, rBj5,
                        SBi, rlSBi, BSi_reg,
                        fSi_theta3, fBi_theta3, fSj_theta3, fBj_theta3,
                        f_Si, f_Bi, f_Sj, f_Bj, f_Bj5);
            }
            // Sj--Bj -> Bi5
            if(Bi5_exists)
            {
                //           Bi5
                //  5'    o -- o===o -- o     3'
                //  ^    /      \        \    |
                //  | P o        \        o P |
                //  |    \        \--.   /    |
                //  |  Si o -- o===o -- o Sj  v
                //  3'       Bi     Bj        5'
                //
                // Bj(5' side) -> Bi5(5' side): cs5_kind

                const auto Bi5_base = parameters_.at(Bi5).base;
                const auto cs_kind  = parameters_.cs5_kind(Bj_base, Bi5_base);
                const auto& rBi5 = sys.position(Bi5);

                energy += calc_force_element(cs_kind, rBj, rBi5,
                        SBj, rlSBj, BSj_reg,
                        fSj_theta3, fBj_theta3, fSi_theta3, fBi_theta3,
                        f_Sj, f_Bj, f_Si, f_Bi, f_Bi5);
            }

            if(Bi3_exists)
            {
                //       Si   Bi   Bj   Sj
                //  5'    o--> o===o <--o     3'
                //  ^    /        /_.'   \    |
                //  | P o        /  tCS   o P |
                //  |    \      /        /    v
                //  3'    o -- o===o -- o     5'
                //           Bi3   Bj5

                const auto Bi3_base = parameters_.at(Bi3).base;
                const auto cs_kind  = parameters_.cs3_kind(Bj_base, Bi3_base);
                const auto& rBi3 = sys.position(Bi3);

                energy += calc_force_element(cs_kind, rBj, rBi3,
                        SBj, rlSBj, BSj_reg,
                        fSj_theta3, fBj_theta3, fSi_theta3, fBi_theta3,
                        f_Sj, f_Bj, f_Si, f_Bi, f_Bi3);

            }
            if(Bj3_exists)
            {
                // ----------------------------------------------------------------
                //
                //           Bi5   Bj3
                //  5'    o--> o===o <--o     3'
                //  ^    /        /_.'   \    |
                //  | P o        /  tCS   o P |
                //  |    \      /        /    v
                //  3'    o -- o===o -- o     5'
                //       Si   Bi   Bj   Sj

                const auto Bj3_base = parameters_.at(Bj3).base;
                const auto cs_kind  = parameters_.cs3_kind(Bi_base, Bj3_base);
                const auto& rBj3 = sys.position(Bj3);

                energy += calc_force_element(cs_kind, rBi, rBj3,
                        SBi, rlSBi, BSi_reg,
                        fSi_theta3, fBi_theta3, fSj_theta3, fBj_theta3,
                        f_Si, f_Bi, f_Sj, f_Bj, f_Bj3);
            }

            sys.force(Si)  += f_Si;
            sys.force(Sj)  += f_Sj;
            sys.force(Bi)  += f_Bi;
            sys.force(Bi3) += f_Bi3;
            sys.force(Bi5) += f_Bi5;
            sys.force(Bj)  += f_Bj;
            sys.force(Bj3) += f_Bj3;
            sys.force(Bj5) += f_Bj5;

            if /*constexpr*/ (NeedVirial)
            {
                sys.virial()  += math::tensor_product(sys.transpose(rSi, rBi), f_Si) +
                                 math::tensor_product(                   rBi , f_Bi) +
                                 math::tensor_product(              rBi + Bij, f_Bj) +
                                 math::tensor_product(sys.transpose(rSj, rBi), f_Sj);
                if(Bi5_exists)
                {
                    sys.virial() += math::tensor_product(sys.transpose(sys.position(Bi5), rBi), f_Bi5);
                }
                if(Bi3_exists)
                {
                    sys.virial() += math::tensor_product(sys.transpose(sys.position(Bi3), rBi), f_Bi3);
                }
                if(Bj5_exists)
                {
                    sys.virial() += math::tensor_product(sys.transpose(sys.position(Bj5), rBi), f_Bj5);
                }
                if(Bj3_exists)
                {
                    sys.virial() += math::tensor_product(sys.transpose(sys.position(Bj3), rBi), f_Bj3);
                }
            }
        }
    }
    return energy;
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{

extern template class ThreeSPN2CrossStackingInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ThreeSPN2CrossStackingInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ThreeSPN2CrossStackingInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ThreeSPN2CrossStackingInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif // MJOLNIR_FORCEFIELD_3SPN_BASE_PAIRING_INTERACTION_HPP
