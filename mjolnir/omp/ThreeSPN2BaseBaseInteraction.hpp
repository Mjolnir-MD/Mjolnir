#ifndef MJOLNIR_OMP_3SPN2_BASE_BASE_INTERACTION_HPP
#define MJOLNIR_OMP_3SPN2_BASE_BASE_INTERACTION_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/UnlimitedGridCellList.hpp>
#include <mjolnir/omp/PeriodicGridCellList.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteraction.hpp>

namespace mjolnir
{

// 3SPN.2 Base-Base non-local interaction (base pairing & cross-stacking).
// This is an implementation of the potential developed in the following paper.
// - D.M.Hinckley, G.S.Freeman, J.K.Whitmer, and J.J.de Pablo (2013) J. Chem. Phys.
//   doi: 10.1063/1.4822042
//
// This interaction is deeply coupled with its potential, so it does not receive
// a potential type as a template argument. It always uses the same potential,
// `ThreeSPN2BaseBaseInteractionPotential`.
//
// It shares CellList between BasePairing and CrossStacking.
// It constructs CellList for BasePairing and re-use the pairs for CrossStacking.
//
template<typename realT, template<typename, typename> class boundaryT>
class ThreeSPN2BaseBaseInteraction<
    OpenMPSimulatorTraits<realT, boundaryT>
    > final : public GlobalInteractionBase<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:

    using traits_type     = OpenMPSimulatorTraits<realT, boundaryT>;
    using base_type       = GlobalInteractionBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using topology_type   = typename base_type::topology_type;
    using boundary_type   = typename base_type::boundary_type;
    using potential_type      = ThreeSPN2BaseBaseInteractionPotential<real_type>;
    using partition_type      = SpatialPartition<traits_type, potential_type>;
    using parameter_list_type = ThreeSPN2BaseBaseInteractionParameterList<traits_type>;

    using base_kind        = parameter_3SPN2::bead_kind;
    using base_pair_kind   = parameter_3SPN2::base_pair_kind;
    using cross_stack_kind = parameter_3SPN2::cross_stack_kind;

  public:

    ThreeSPN2BaseBaseInteraction(parameter_list_type&& para, partition_type&& part)
        : parameters_(std::move(para)), partition_(std::move(part))
    {}
    ~ThreeSPN2BaseBaseInteraction() override {}

    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->parameters_.initialize(sys, topol);
        this->partition_ .initialize(sys, this->parameters_);
    }

    void update(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->parameters_.update(sys, topol);
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

    void      calc_force (system_type& sys)       const noexcept override
    {
        constexpr auto pi        = math::constants<real_type>::pi();
        constexpr auto two_pi    = math::constants<real_type>::two_pi();
        constexpr auto tolerance = math::abs_tolerance<real_type>();

        const auto leading_participants = this->parameters_.leading_participants();
#pragma omp parallel for
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
        {
            const std::size_t thread_id = omp_get_thread_num();

            const auto   Bi = leading_participants[idx];
            const auto& rBi = sys.position(Bi);
            for(const auto& ptnr : this->partition_.partners(Bi))
            {
                const auto  Bj  = ptnr.index;
                const auto& pot = ptnr.potential();
                const auto& rBj = sys.position(Bj);

                const auto Bij = sys.adjust_direction(rBi, rBj); // Bi -> Bj
                const auto lBij_sq = math::length_sq(Bij);
                if(lBij_sq > parameters_.cutoff_sq())
                {
                    continue;
                }

                // ================================================================
                // base pairing
                //
                //  Si o         o Sj
                //      \-.   ,-/
                //    Bi o =(= o Bj
                //
                // U_rep(rij) + 1/2(1+cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)

                const auto rlBij = math::rsqrt(lBij_sq); // 1 / |Bij|
                const auto lBij  = lBij_sq * rlBij;      // |Bij|

                const auto Bij_reg =  rlBij * Bij;
                const auto Bji_reg = -rlBij * Bij;

                const auto bp_kind    = pot.bp_kind;
                const auto epsilon_bp = parameters_.epsilon(bp_kind);
                const auto alpha_bp   = parameters_.alpha(bp_kind);
                const auto r0_bp      = parameters_.r0(bp_kind);

                // ----------------------------------------------------------------
                // calculate the the repulsive part, which does not depend on angle.
                //
                // dU_rep = 2 a e exp(-a(r-r0)) (1-exp(-a(r-r0))) ... r  <  r0
                //        = 0                                     ... r0 <= r
                //
                {
                    const auto dU_rep = pot.dU_rep(epsilon_bp, alpha_bp, lBij, r0_bp);
                    if(dU_rep != real_type(0))
                    {
                        // remember that F = -dU.
                        const auto F = -dU_rep * Bji_reg;
                        sys.force_thread(thread_id, Bi) += F;
                        sys.force_thread(thread_id, Bj) -= F;

                        // Bij = Bi -> Bj = rBj - rBi
                        sys.virial_thread(thread_id) += math::tensor_product(Bij, -F);
                    }
                }

                // ----------------------------------------------------------------
                // calc theta1 and 2 to calculate the attractive part,
                //  = 1/2(1+cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)
                //
                //   theta1   theta2
                //       |     |
                //  Si o v     v o Sj
                //      \-.   ,-/
                //    Bi o =(= o Bj

                const auto   Si = pot.Si;
                const auto   Sj = pot.Sj;
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

                const auto dot_SBiBj  = -math::dot_product(SBi, Bij);
                const auto dot_SBjBi  =  math::dot_product(SBj, Bij);
                const auto cos_theta1 = dot_SBiBj * rlSBi * rlBij;
                const auto cos_theta2 = dot_SBjBi * rlSBj * rlBij;
                const auto theta1 = std::acos(math::clamp<real_type>(cos_theta1, -1, 1));
                const auto theta2 = std::acos(math::clamp<real_type>(cos_theta2, -1, 1));

                // ----------------------------------------------------------------
                // calc angle-dependent terms and advance if both are nonzero
                //
                // 1/2(1+cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)

                const auto K_BP         = parameters_.K_BP();
                const auto pi_over_K_BP = parameters_.pi_over_K_BP();
                const auto theta1_0     = parameters_.theta1_0(bp_kind);
                const auto theta2_0     = parameters_.theta2_0(bp_kind);

                const auto f1 = pot.f(K_BP, pi_over_K_BP, theta1, theta1_0);
                const auto f2 = pot.f(K_BP, pi_over_K_BP, theta2, theta2_0);

                if(f1 != real_type(0.0) && f2 != real_type(0.0))
                {
                    // calculate dihedral, phi
                    //
                    //  Si o         o Sj
                    //      \       /
                    //    Bi o =(= o Bj
                    //         phi

                    const auto df1 = pot.df(K_BP, pi_over_K_BP, theta1, theta1_0);
                    const auto df2 = pot.df(K_BP, pi_over_K_BP, theta2, theta2_0);

                    const auto m = math::cross_product(-SBi, Bij);
                    const auto n = math::cross_product( Bij, SBj);
                    const auto m_lsq = math::length_sq(m);
                    const auto n_lsq = math::length_sq(n);

                    const auto dot_phi = math::dot_product(m, n) *
                                         math::rsqrt(m_lsq * n_lsq);
                    const auto cos_phi = math::clamp<real_type>(dot_phi, -1, 1);

                    const auto phi = std::copysign(std::acos(cos_phi),
                                                   -math::dot_product(SBi, n));

                    auto dphi = phi - this->parameters_.phi_0(bp_kind);
                    if(dphi < -pi) {dphi += two_pi;}
                    if(pi <= dphi) {dphi -= two_pi;}
                    const auto cos_dphi = std::cos(dphi);
                    const auto sin_dphi = std::sin(dphi);

                    // ------------------------------------------------------------
                    // calculate attractive force
                    //
                    // d/dr [1/2 (1 + cos(dphi)) f(dtheta1) f(dtheta2) U_attr(Bij)]
                    // = ( -sin(dphi))/2 f(dtheta1) f(dtheta2) U_attr(Bij) dphi/dr
                    // + (1+cos(dphi))/2 df/dtheta1 f(dtheta2) U_attr(Bij) dtheta1/dr
                    // + (1+cos(dphi))/2 f(dtheta1) df/dtheta2 U_attr(Bij) dtheta2/dr
                    // + (1+cos(dphi))/2 f(dtheta1) f(dtheta2) dU_attr/dr  dBij/dr

                    if(cos_dphi != real_type(-1.0))
                    {
                        auto f_Si = math::make_coordinate<coordinate_type>(0,0,0);
                        auto f_Sj = math::make_coordinate<coordinate_type>(0,0,0);
                        auto f_Bi = math::make_coordinate<coordinate_type>(0,0,0);
                        auto f_Bj = math::make_coordinate<coordinate_type>(0,0,0);

                        const auto U_dU_attr = pot.U_dU_attr(epsilon_bp, alpha_bp, lBij, r0_bp);

                        // --------------------------------------------------------
                        // calc dihedral term
                        // -sin(dphi)/2 f(dtheta1) f(dtheta2) U_attr(Bij) dphi/dr
                        if(sin_dphi != real_type(0.0))
                        {
                            // remember that F = -dU. Here `coef` = -dU.
                            const auto coef = real_type(0.5) * sin_dphi *
                                              f1 * f2 * U_dU_attr.first;

                            const auto rlBij_sq = rlBij * rlBij; // 1 / |Bij|^2
                            const auto fSi = ( coef * lBij / m_lsq) * m;
                            const auto fSj = (-coef * lBij / n_lsq) * n;

                            const auto coef_Bi = dot_SBiBj * rlBij_sq;
                            const auto coef_Bj = dot_SBjBi * rlBij_sq;

                            f_Si += fSi;
                            f_Bi += (coef_Bi - real_type(1.0)) * fSi - coef_Bj * fSj;
                            f_Bj += (coef_Bj - real_type(1.0)) * fSj - coef_Bi * fSi;
                            f_Sj += fSj;
                        }

                        const auto dihd_term = real_type(0.5) * (real_type(1.0) + cos_dphi);

                        // --------------------------------------------------------
                        // calc theta1 term
                        // (1+cos(dphi))/2 df/dtheta1 f(dtheta2) U_attr(Bij) dtheta1/dr
                        if(df1 != real_type(0.0))
                        {
                            // remember that F = -dU. Here `coef` = -dU.
                            const auto coef = -dihd_term * df1 * f2 * U_dU_attr.first;

                            const auto sin_theta1 = std::sin(theta1);
                            const auto coef_rsin  = (sin_theta1 > tolerance) ?
                                       (coef / sin_theta1) : (coef / tolerance);

                            const auto fSi = (coef_rsin * rlSBi) *
                                             (cos_theta1 * BSi_reg - Bij_reg);
                            const auto fBj = (coef_rsin * rlBij) *
                                             (cos_theta1 * Bij_reg - BSi_reg);
                            f_Si += fSi;
                            f_Bi -= (fSi + fBj);
                            f_Bj += fBj;
                        }
                        // --------------------------------------------------------
                        // calc theta2 term
                        // (1+cos(dphi))/2 f(dtheta1) df/dtheta2 U_attr(Bij) dtheta2/dr
                        if(df2 != real_type(0.0))
                        {
                            // remember that F = -dU. Here `coef` = -dU.
                            const auto coef = -dihd_term * f1 * df2 * U_dU_attr.first;

                            const auto sin_theta2 = std::sin(theta2);
                            const auto coef_rsin  = (sin_theta2 > tolerance) ?
                                       (coef / sin_theta2) : (coef / tolerance);

                            const auto fBi = (coef_rsin * rlBij) *
                                             (cos_theta2 * Bji_reg - BSj_reg);
                            const auto fSj = (coef_rsin * rlSBj) *
                                             (cos_theta2 * BSj_reg - Bji_reg);
                            f_Bi += fBi;
                            f_Bj -= (fBi + fSj);
                            f_Sj += fSj;
                        }
                        // --------------------------------------------------------
                        // calc distance
                        // + 1/2 (1+cos(dphi)) f(dtheta1) f(dtheta2) dU_attr/dr  dBij/dr
                        if(U_dU_attr.second != real_type(0.0))
                        {
                            // remember that F = -dU. Here `coef` = -dU.
                            const auto coef = -dihd_term * f1 * f2 * U_dU_attr.second;
                            f_Bi += coef * Bji_reg;
                            f_Bj += coef * Bij_reg;
                        }

                        sys.force_thread(thread_id, Si) += f_Si;
                        sys.force_thread(thread_id, Bi) += f_Bi;
                        sys.force_thread(thread_id, Bj) += f_Bj;
                        sys.force_thread(thread_id, Sj) += f_Sj;

                        sys.virial_thread(thread_id)  +=
                            math::tensor_product(sys.transpose(rSi, rBi), f_Si) +
                            math::tensor_product(                   rBi , f_Bi) +
                            math::tensor_product(             rBi + Bij , f_Bj) +
                            math::tensor_product(sys.transpose(rSj, rBi), f_Sj) ;
                    }
                }

                // ================================================================
                // cross stacking
                // f(theta_3) f(theta_CS) U_attr(epsilon, alpha, rij)
                //
                //       Si   Bi   Bj   Sj
                //  5'    o -- o===o -- o     3'
                //  ^    /      \ /      \    |
                //  | P o        X        o P |
                //  |    \      / \      /    v
                //  3'    o -- o===o -- o     5'
                //       Bj_next   Bj_next
                //
                // d/dr Vcs =
                //    df/dtheta3 f(theta_CS)  U_attr(eps, alp, rij) dtheta_3  /dr
                //  + f(theta_3) df/dtheta_CS U_attr(eps, alp, rij) dtheta_CS /dr
                //  + f(theta_3) f(theta_CS)  dU_attr/drij          drij/dr
                //

                const auto Bi_next        = pot.Bi_next;
                const auto Bj_next        = pot.Bj_next;
                const bool Bi_next_exists = (Bi_next != potential_type::invalid());
                const bool Bj_next_exists = (Bj_next != potential_type::invalid());

                if(!Bi_next_exists && !Bj_next_exists)
                {
                    continue; // if both interacting pair do not exist, do nothing.
                }

                // ----------------------------------------------------------------
                // calc common part, theta3 and dtheta3/dr.

                const auto cos_theta3 = math::dot_product(BSi_reg, BSj_reg);
                const auto theta3     = std::acos(math::clamp<real_type>(cos_theta3, -1, 1));
                const auto theta3_0   = parameters_.theta3_0(bp_kind);
                const auto f3         = pot.f(K_BP, pi_over_K_BP, theta3, theta3_0);
                if(f3 == real_type(0))
                {
                    // f(theta) == 0 means df(theta) is also zero.
                    // so here, both cross-stacking becomes zero. skip them.
                    continue;
                }
                const auto df3 = pot.df(K_BP, pi_over_K_BP, theta3, theta3_0);

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

                const auto K_CS            = parameters_.K_CS();
                const auto pi_over_K_CS    = parameters_.pi_over_K_CS();

                auto f_Si      = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Sj      = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Bi      = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Bj      = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Bi_next = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Bj_next = math::make_coordinate<coordinate_type>(0,0,0);

                // adjacent of Base j exists. its not the edge of the strand.
                if(Bj_next_exists)
                {
                    // ------------------------------------------------------------
                    // 5' cross stacking (the case if `Bi` is in sense strand)
                    //
                    //       Si   Bi   Bj   Sj
                    //  5'    o--> o===o <--o     3'
                    //  ^    /   `--\        \    |
                    //  | P o   tCS  \        o P |
                    //  |    \        \      /    v
                    //  3'    o -- o===o -- o     5'
                    //           Bi3   Bj5
                    //
                    // Hereafter, we use the notation illustrated above (excepting
                    // the index `Bj_next` that corresponds to Bj5).

                    const auto cs_kind = pot.cs_i_kind;

                    const auto epsilon_cs = parameters_.epsilon(cs_kind);
                    const auto alpha_cs   = parameters_.alpha(cs_kind);
                    const auto r0_cs      = parameters_.r0(cs_kind);

                    const auto& rBj5     = sys.position(Bj_next);
                    const auto  Bj5i     = sys.adjust_direction(rBj5, rBi);
                    const auto  lBj5i_sq = math::length_sq(Bj5i); // |Bj5i|^2
                    const auto  rlBj5i   = math::rsqrt(lBj5i_sq);  // 1 / |Bj5i|

                    const auto dot_theta_CS = math::dot_product(SBi, Bj5i);
                    const auto cos_theta_CS = dot_theta_CS * rlSBi * rlBj5i;
                    const auto theta_CS     =
                        std::acos(math::clamp<real_type>(cos_theta_CS, -1, 1));
                    const auto theta_CS_0   = parameters_.thetaCS_0(cs_kind);

                    const auto fCS  = pot.f(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                    // if f == 0, df is also zero. if fCS == 0, no force there
                    if(fCS != real_type(0))
                    {
                        const auto dfCS  = pot.df(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);

                        const auto lBj5i = lBj5i_sq * rlBj5i;
                        const auto U_dU_attr = pot.U_dU_attr(epsilon_cs, alpha_cs, lBj5i, r0_cs);
                        const auto Bj5i_reg  = rlBj5i * Bj5i;

                        // --------------------------------------------------------
                        // df/dtheta3 f(theta_CS)  U_attr(eps, alp, rij) dtheta_3 /dr
                        if(df3 != real_type(0))
                        {
                            const auto coef = -df3 * fCS * U_dU_attr.first;
                            f_Si += coef * fSi_theta3;
                            f_Sj += coef * fSj_theta3;
                            f_Bi += coef * fBi_theta3;
                            f_Bj += coef * fBj_theta3;
                        }
                        // --------------------------------------------------------
                        // f(theta_3) df/dtheta_CS U_attr(eps, alp, rij) dtheta_CS/dr
                        if(dfCS != real_type(0))
                        {
                            const auto coef         = -f3 * dfCS * U_dU_attr.first;
                            const auto sin_theta_CS = std::sin(theta_CS);
                            const auto coef_rsin    = (sin_theta_CS > tolerance) ?
                                       (coef / sin_theta_CS) : (coef / tolerance);

                            const auto fSi  =  coef_rsin * rlSBi  * (cos_theta_CS * BSi_reg + Bj5i_reg);
                            const auto fBj5 = -coef_rsin * rlBj5i * (cos_theta_CS * Bj5i_reg + BSi_reg);

                            f_Si      += fSi;
                            f_Bi      -= (fSi + fBj5);
                            f_Bj_next += fBj5;
                        }
                        // --------------------------------------------------------
                        // f(theta_3) f(theta_CS)  dU_attr/drij          drij/dr
                        if(U_dU_attr.second != real_type(0.0))
                        {
                            const auto coef = -f3 * fCS * U_dU_attr.second;
                            f_Bi      += coef * Bj5i_reg;
                            f_Bj_next -= coef * Bj5i_reg;
                        }
                    }
                }
                if(Bi_next_exists)
                {
                    // ------------------------------------------------------------
                    // 3' cross stacking (the case `Bi` is in a sense strand)
                    //
                    //       Si   Bi   Bj   Sj
                    //  5'    o--> o===o <--o     3'
                    //  ^    /        /--'   \    |
                    //  | P o        /  tCS   o P |
                    //  |    \      /        /    v
                    //  3'    o -- o===o -- o     5'
                    //           Bi3   Bj5
                    //
                    // Hereafter, we use the notation illustrated above (excepting
                    // the index `Bi_next` that corresponds to Bi3).

                    const auto cs_kind = pot.cs_j_kind;

                    const auto epsilon_cs = parameters_.epsilon(cs_kind);
                    const auto alpha_cs   = parameters_.alpha(cs_kind);
                    const auto r0_cs      = parameters_.r0(cs_kind);

                    const auto& rBi3     = sys.position(Bi_next);
                    const auto  Bi3j     = sys.adjust_direction(rBi3, rBj);
                    const auto  lBi3j_sq = math::length_sq(Bi3j); // |Bi3j|^2
                    const auto  rlBi3j   = math::rsqrt(lBi3j_sq);  // 1 / |Bi3j|

                    const auto dot_theta_CS = math::dot_product(SBj, Bi3j);
                    const auto cos_theta_CS = dot_theta_CS * rlSBj * rlBi3j;
                    const auto theta_CS     =
                        std::acos(math::clamp<real_type>(cos_theta_CS, -1, 1));
                    const auto theta_CS_0   = parameters_.thetaCS_0(cs_kind);

                    const auto fCS = pot.f(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                    // if f == 0, df is also zero. if fCS == 0, no force there
                    if(fCS != real_type(0))
                    {
                        const auto dfCS  = pot.df(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                        const auto lBi3j = lBi3j_sq * rlBi3j;
                        const auto U_dU_attr = pot.U_dU_attr(epsilon_cs, alpha_cs, lBi3j, r0_cs);
                        const auto Bi3j_reg  = rlBi3j * Bi3j;

                        // --------------------------------------------------------
                        // df/dtheta3 f(theta_CS)  U_attr(eps, alp, rij) dtheta_3 /dr
                        if(df3 != real_type(0))
                        {
                            const auto coef = -df3 * fCS * U_dU_attr.first;
                            f_Si += coef * fSi_theta3;
                            f_Sj += coef * fSj_theta3;
                            f_Bi += coef * fBi_theta3;
                            f_Bj += coef * fBj_theta3;
                        }
                        // --------------------------------------------------------
                        // f(theta_3) df/dtheta_CS U_attr(eps, alp, rij) dtheta_CS/dr
                        if(dfCS != real_type(0))
                        {
                            const auto coef         = -f3 * dfCS * U_dU_attr.first;
                            const auto sin_theta_CS = std::sin(theta_CS);
                            const auto coef_rsin    = (sin_theta_CS > tolerance) ?
                                       (coef / sin_theta_CS) : (coef / tolerance);

                            const auto fSj  =  coef_rsin * rlSBj  * (cos_theta_CS * BSj_reg + Bi3j_reg);
                            const auto fBi3 = -coef_rsin * rlBi3j * (cos_theta_CS * Bi3j_reg + BSj_reg);
                            f_Sj      += fSj;
                            f_Bj      -= (fSj + fBi3);
                            f_Bi_next += fBi3;
                        }
                        // --------------------------------------------------------
                        // f(theta_3) f(theta_CS)  dU_attr/drij          drij/dr
                        if(U_dU_attr.second != real_type(0.0))
                        {
                            const auto coef = -f3 * fCS * U_dU_attr.second;
                            f_Bj      += coef * Bi3j_reg;
                            f_Bi_next -= coef * Bi3j_reg;
                        }
                    }
                }

                sys.force_thread(thread_id, Si)      += f_Si     ;
                sys.force_thread(thread_id, Sj)      += f_Sj     ;
                sys.force_thread(thread_id, Bi)      += f_Bi     ;
                sys.force_thread(thread_id, Bj)      += f_Bj     ;
                sys.force_thread(thread_id, Bi_next) += f_Bi_next;
                sys.force_thread(thread_id, Bj_next) += f_Bj_next;

                sys.virial_thread(thread_id)  +=
                    math::tensor_product(sys.transpose(rSi, rBi), f_Si) +
                    math::tensor_product(                   rBi , f_Bi) +
                    math::tensor_product(              rBi + Bij, f_Bj) +
                    math::tensor_product(sys.transpose(rSj, rBi), f_Sj) +
                    math::tensor_product(sys.transpose(sys.position(Bi_next), rBi), f_Bi_next) +
                    math::tensor_product(sys.transpose(sys.position(Bj_next), rBi), f_Bj_next) ;
            }
        }
        return;
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        constexpr auto pi        = math::constants<real_type>::pi();
        constexpr auto two_pi    = math::constants<real_type>::two_pi();

        real_type E = 0.0;

        const auto leading_participants = this->parameters_.leading_participants();
#pragma omp parallel for reduction(+:E)
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
        {
            const auto   Bi = leading_participants[idx];
            const auto& rBi = sys.position(Bi);

            for(const auto& ptnr : this->partition_.partners(Bi))
            {
                const auto  Bj  = ptnr.index;
                const auto& pot = ptnr.potential();
                const auto& rBj = sys.position(Bj);

                const auto  bp_kind = pot.bp_kind;

                const auto Bij = sys.adjust_direction(rBi, rBj); // Bi -> Bj

                const auto lBij_sq = math::length_sq(Bij);
                if(lBij_sq > parameters_.cutoff_sq())
                {
                    continue;
                }

                // ----------------------------------------------------------------
                // base pairing
                //
                //  Si o         o Sj
                //      \-.   ,-/
                //    Bi o =(= o Bj
                //
                // U_rep(rij) + 1/2(1+cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)

                const auto rlBij = math::rsqrt(lBij_sq); // 1 / |Bij|
                const auto lBij  = lBij_sq * rlBij;      // |Bij|

                // ----------------------------------------------------------------
                // calc U_rep. it always exists because it does not depend on angle.
                //
                // U_rep = e_ij(1 - exp(-a_ij(rij - r0_ij)))^2 ... r < r0
                //         0                                   ... otherwise

                const real_type epsilon_bp = parameters_.epsilon(bp_kind);
                const real_type alpha_bp   = parameters_.alpha(bp_kind);
                const real_type r0_bp      = parameters_.r0(bp_kind);

                E += pot.U_rep(epsilon_bp, alpha_bp, lBij, r0_bp);

                // ----------------------------------------------------------------
                // calc theta1 and 2
                //
                //   theta1  theta2
                //       |    |
                //  Si o v    v o Sj
                //      \-.  ,-/
                //    Bi o == o Bj

                const auto   Si = pot.Si;
                const auto   Sj = pot.Sj;
                const auto& rSi = sys.position(Si);
                const auto& rSj = sys.position(Sj);

                const auto SBi = sys.adjust_direction(rSi, rBi); // Si -> Bi
                const auto SBj = sys.adjust_direction(rSj, rBj); // Sj -> Bj

                const auto lSBi_sq = math::length_sq(SBi); // |SBi|^2
                const auto lSBj_sq = math::length_sq(SBj); // |SBj|^2

                const auto rlSBi = math::rsqrt(lSBi_sq); // 1 / |SBi|
                const auto rlSBj = math::rsqrt(lSBj_sq); // 1 / |SBj|

                const auto dot_SBiBj = -math::dot_product(SBi, Bij);
                const auto dot_SBjBi =  math::dot_product(SBj, Bij);

                const auto cos_theta1 = dot_SBiBj * rlSBi * rlBij;
                const auto cos_theta2 = dot_SBjBi * rlSBj * rlBij;

                const auto theta1 = std::acos(math::clamp<real_type>(cos_theta1, -1, 1));
                const auto theta2 = std::acos(math::clamp<real_type>(cos_theta2, -1, 1));

                // ----------------------------------------------------------------
                // The second term of base-pairing
                //  = 1/2(1+cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)
                //
                // f(theta) = 1                             ... abs(dtheta) < pi/2K
                //            1 - cos^2(K (theta - theta0)) ... pi/2K < abs(dtheta) < pi/K
                //            0                             ... pi/K  < abs(dtheta)

                const auto K_BP         = parameters_.K_BP();
                const auto pi_over_K_BP = parameters_.pi_over_K_BP();
                const auto theta1_0     = parameters_.theta1_0(bp_kind);
                const auto theta2_0     = parameters_.theta2_0(bp_kind);

                const auto f1 = pot.f(K_BP, pi_over_K_BP, theta1, theta1_0);
                const auto f2 = pot.f(K_BP, pi_over_K_BP, theta2, theta2_0);

                if(f1 != real_type(0.0) && f2 != real_type(0.0)) // [[likely]]
                {
                    // if both fs are non-zero, the attractive part of base-pairing
                    // has a non-zero value. calculate dihedral and cos(dphi).

                    // ------------------------------------------------------------
                    //  Si o         o Si
                    //      \       /
                    //    Bi o =(= o Bj
                    //         phi

                    const auto Bij_reg = Bij * rlBij;
                    const auto R = math::dot_product(SBi, Bij_reg) * Bij_reg - SBi;
                    const auto S = math::dot_product(SBj, Bij_reg) * Bij_reg - SBj;

                    const auto dot_RS  = math::dot_product(R, S) *
                            math::rsqrt(math::length_sq(R) * math::length_sq(S));
                    const auto cos_phi = math::clamp<real_type>(dot_RS, -1, 1);

                    const auto n    =  math::cross_product(Bij, SBj);
                    const auto sign = -math::dot_product(SBi, n);
                    const auto phi  = std::copysign(std::acos(cos_phi), sign);

                    auto dphi = phi - parameters_.phi_0(bp_kind);
                    if     (pi   < dphi) {dphi -= two_pi;}
                    else if(dphi <  -pi) {dphi += two_pi;}

                    const auto cos_dphi = std::cos(dphi);

                    // ------------------------------------------------------------
                    // U_attr = -e_ij                                       .. r < r0
                    //          -e_ij + e_ij(1 - exp(-a_ij(rij - r0_ij)))^2 .. r0 < r
                    //
                    // XXX: Note that std::exp would be calculated only once because
                    //      the condition that requires exp differs.
                    const auto U_attr = pot.U_attr(epsilon_bp, alpha_bp, lBij, r0_bp);

                    // ------------------------------------------------------------
                    // The second term of base-pairing
                    //  = 1/2(1 + cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)
                    E += 0.5 * (1 + cos_dphi) * f1 * f2 * U_attr;
                }

                // ================================================================
                // cross stacking
                // f(theta_3) f(theta_CS) U_attr(epsilon, alpha, rij)
                //
                //       Si   Bi   Bj   Sj
                //  5'    o -- o===o -- o     3'
                //  ^    /      \ /      \    |
                //  | P o        x        o P |
                //  |    \      / \      /    v
                //  3'    o -- o===o -- o     5'
                //       Bj_next   Bj_next

                const auto Bi_next        = pot.Bi_next;
                const auto Bj_next        = pot.Bj_next;
                const bool Bi_next_exists = (Bi_next != potential_type::invalid());
                const bool Bj_next_exists = (Bj_next != potential_type::invalid());

                if(!Bi_next_exists && !Bj_next_exists)
                {
                    continue; // if both interacting pair do not exist, do nothing.
                }

                const auto dot_SBiSBj = math::dot_product(SBi, SBj);
                const auto cos_theta3 = dot_SBiSBj * rlSBi * rlSBj;
                const auto theta3     = std::acos(math::clamp<real_type>(cos_theta3, -1, 1));
                const auto theta3_0   = parameters_.theta3_0(bp_kind);
                const auto f3         = pot.f(K_BP, pi_over_K_BP, theta3, theta3_0);
                if(f3 == real_type(0))
                {
                    continue; // both cross-stacking becomes zero. skip them.
                }

                const auto K_CS         = parameters_.K_CS();
                const auto pi_over_K_CS = parameters_.pi_over_K_CS();

                // adjacent of Base Bj exists. its not the edge of the strand.
                if(Bj_next_exists)
                {
                    // ----------------------------------------------------------------
                    // cross stacking
                    // f(theta_3) f(theta_CS) U_attr(epsilon, alpha, rij)
                    //
                    //       Si   Bi   Bj   Sj
                    //  5'    o--> o===o <--o     3'
                    //  ^    /   `--\        \    |
                    //  | P o   tCS  \        o P |
                    //  |    \        \      /    v
                    //  3'    o -- o===o -- o     5'
                    //           Bi3   Bj5
                    //
                    // Hereafter, we use the notation illustrated above (excepting
                    // the index `Bj_next` that corresponds to Bj5).

                    const auto cs_kind    = pot.cs_i_kind;
                    const auto epsilon_cs = parameters_.epsilon(cs_kind);
                    const auto alpha_cs   = parameters_.alpha(cs_kind);
                    const auto r0_cs      = parameters_.r0(cs_kind);

                    const auto& rBj5 = sys.position(Bj_next);

                    const auto Bj5i     = sys.adjust_direction(rBj5, rBi);
                    const auto lBj5i_sq = math::length_sq(Bj5i);
                    const auto rlBj5i   = math::rsqrt(lBj5i_sq);

                    const auto dot_theta_CS = math::dot_product(SBi, Bj5i);
                    const auto cos_theta_CS = dot_theta_CS * rlSBi * rlBj5i;
                    const auto theta_CS     = std::acos(
                            math::clamp<real_type>(cos_theta_CS, -1, 1));
                    const auto theta_CS_0   = parameters_.thetaCS_0(cs_kind);

                    const auto fCS = pot.f(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                    if(fCS != real_type(0))
                    {
                        const auto lBj5i  = lBj5i_sq * rlBj5i;
                        const auto U_attr = pot.U_attr(epsilon_cs, alpha_cs, lBj5i, r0_cs);
                        E += f3 * fCS * U_attr;
                    }
                }
                if(Bi_next_exists)
                {
                    // ----------------------------------------------------------------
                    // cross stacking
                    // f(theta_3) f(theta_CS) U_attr(epsilon, alpha, rij)
                    //
                    //       Si   Bi   Bj   Sj
                    //  5'    o--> o===o <--o     3'
                    //  ^    /        /--'   \    |
                    //  | P o        /  tCS   o P |
                    //  |    \      /        /    v
                    //  3'    o -- o===o -- o     5'
                    //           Bi3   Bj5
                    //
                    // Hereafter, we use the notation illustrated above (excepting
                    // the index `Bi_next` that corresponds to Bi3).

                    const auto cs_kind    = pot.cs_j_kind;
                    const auto epsilon_cs = parameters_.epsilon(cs_kind);
                    const auto alpha_cs   = parameters_.alpha(cs_kind);
                    const auto r0_cs      = parameters_.r0(cs_kind);

                    const auto& rBi3    = sys.position(Bi_next);

                    const auto Bi3j     = sys.adjust_direction(rBi3, rBj);
                    const auto lBi3j_sq = math::length_sq(Bi3j);
                    const auto rlBi3j   = math::rsqrt(lBi3j_sq);

                    const auto dot_theta_CS = math::dot_product(SBj, Bi3j);
                    const auto cos_theta_CS = dot_theta_CS * rlSBj * rlBi3j;
                    const auto theta_CS     = std::acos(
                            math::clamp<real_type>(cos_theta_CS, -1, 1));
                    const auto theta_CS_0   = parameters_.thetaCS_0(cs_kind);

                    const auto fCS = pot.f(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                    if(fCS != real_type(0))
                    {
                        const auto lBi3j  = lBi3j_sq * rlBi3j;
                        const auto U_attr = pot.U_attr(epsilon_cs, alpha_cs, lBi3j, r0_cs);
                        E += f3 * fCS * U_attr;
                    }
                }
            }
        }
        return E;
    }

    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        constexpr auto pi        = math::constants<real_type>::pi();
        constexpr auto two_pi    = math::constants<real_type>::two_pi();
        constexpr auto tolerance = math::abs_tolerance<real_type>();

        const auto leading_participants = this->parameters_.leading_participants();
        real_type energy = 0;
#pragma omp parallel for reduction(+:energy)
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
        {
            const std::size_t thread_id = omp_get_thread_num();

            const auto   Bi = leading_participants[idx];
            const auto& rBi = sys.position(Bi);
            for(const auto& ptnr : this->partition_.partners(Bi))
            {
                const auto  Bj  = ptnr.index;
                const auto& pot = ptnr.potential();
                const auto& rBj = sys.position(Bj);

                const auto Bij = sys.adjust_direction(rBi, rBj); // Bi -> Bj
                const auto lBij_sq = math::length_sq(Bij);
                if(lBij_sq > parameters_.cutoff_sq())
                {
                    continue;
                }

                // ================================================================
                // base pairing
                //
                //  Si o         o Sj
                //      \-.   ,-/
                //    Bi o =(= o Bj
                //
                // U_rep(rij) + 1/2(1+cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)

                const auto rlBij = math::rsqrt(lBij_sq); // 1 / |Bij|
                const auto lBij  = lBij_sq * rlBij;      // |Bij|

                const auto Bij_reg =  rlBij * Bij;
                const auto Bji_reg = -rlBij * Bij;

                const auto bp_kind    = pot.bp_kind;
                const auto epsilon_bp = parameters_.epsilon(bp_kind);
                const auto alpha_bp   = parameters_.alpha(bp_kind);
                const auto r0_bp      = parameters_.r0(bp_kind);

                // ----------------------------------------------------------------
                // calculate the the repulsive part, which does not depend on angle.
                //
                // dU_rep = 2 a e exp(-a(r-r0)) (1-exp(-a(r-r0))) ... r  <  r0
                //        = 0                                     ... r0 <= r
                //
                {
                    energy += pot.U_rep(epsilon_bp, alpha_bp, lBij, r0_bp);
                    const auto dU_rep = pot.dU_rep(epsilon_bp, alpha_bp, lBij, r0_bp);
                    if(dU_rep != real_type(0))
                    {
                        // remember that F = -dU.
                        const auto F = -dU_rep * Bji_reg;
                        sys.force_thread(thread_id, Bi) += F;
                        sys.force_thread(thread_id, Bj) -= F;

                        sys.virial_thread(thread_id) += math::tensor_product(Bij, -F); // (Bj - Bi) * Fj
                    }
                }

                // ----------------------------------------------------------------
                // calc theta1 and 2 to calculate the attractive part,
                //  = 1/2(1+cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)
                //
                //   theta1   theta2
                //       |     |
                //  Si o v     v o Sj
                //      \-.   ,-/
                //    Bi o =(= o Bj

                const auto   Si = pot.Si;
                const auto   Sj = pot.Sj;
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

                const auto dot_SBiBj  = -math::dot_product(SBi, Bij);
                const auto dot_SBjBi  =  math::dot_product(SBj, Bij);
                const auto cos_theta1 = dot_SBiBj * rlSBi * rlBij;
                const auto cos_theta2 = dot_SBjBi * rlSBj * rlBij;
                const auto theta1 = std::acos(math::clamp<real_type>(cos_theta1, -1, 1));
                const auto theta2 = std::acos(math::clamp<real_type>(cos_theta2, -1, 1));

                // ----------------------------------------------------------------
                // calc angle-dependent terms and advance if both are nonzero
                //
                // 1/2(1+cos(dphi)) f(dtheta1) f(dtheta2) U_attr(rij)

                const auto K_BP         = parameters_.K_BP();
                const auto pi_over_K_BP = parameters_.pi_over_K_BP();
                const auto theta1_0     = parameters_.theta1_0(bp_kind);
                const auto theta2_0     = parameters_.theta2_0(bp_kind);

                const auto f1 = pot.f(K_BP, pi_over_K_BP, theta1, theta1_0);
                const auto f2 = pot.f(K_BP, pi_over_K_BP, theta2, theta2_0);

                if(f1 != real_type(0.0) && f2 != real_type(0.0))
                {
                    // calculate dihedral, phi
                    //
                    //  Si o         o Sj
                    //      \       /
                    //    Bi o =(= o Bj
                    //         phi

                    const auto df1 = pot.df(K_BP, pi_over_K_BP, theta1, theta1_0);
                    const auto df2 = pot.df(K_BP, pi_over_K_BP, theta2, theta2_0);

                    const auto m = math::cross_product(-SBi, Bij);
                    const auto n = math::cross_product( Bij, SBj);
                    const auto m_lsq = math::length_sq(m);
                    const auto n_lsq = math::length_sq(n);

                    const auto dot_phi = math::dot_product(m, n) *
                                         math::rsqrt(m_lsq * n_lsq);
                    const auto cos_phi = math::clamp<real_type>(dot_phi, -1, 1);

                    const auto phi = std::copysign(std::acos(cos_phi),
                                                   -math::dot_product(SBi, n));

                    auto dphi = phi - this->parameters_.phi_0(bp_kind);
                    if(dphi < -pi) {dphi += two_pi;}
                    if(pi <= dphi) {dphi -= two_pi;}
                    const auto cos_dphi = std::cos(dphi);
                    const auto sin_dphi = std::sin(dphi);

                    // ------------------------------------------------------------
                    // calculate attractive force
                    //
                    // d/dr [1/2 (1 + cos(dphi)) f(dtheta1) f(dtheta2) U_attr(Bij)]
                    // = ( -sin(dphi))/2 f(dtheta1) f(dtheta2) U_attr(Bij) dphi/dr
                    // + (1+cos(dphi))/2 df/dtheta1 f(dtheta2) U_attr(Bij) dtheta1/dr
                    // + (1+cos(dphi))/2 f(dtheta1) df/dtheta2 U_attr(Bij) dtheta2/dr
                    // + (1+cos(dphi))/2 f(dtheta1) f(dtheta2) dU_attr/dr  dBij/dr

                    if(cos_dphi != real_type(-1.0))
                    {
                        auto f_Si = math::make_coordinate<coordinate_type>(0,0,0);
                        auto f_Sj = math::make_coordinate<coordinate_type>(0,0,0);
                        auto f_Bi = math::make_coordinate<coordinate_type>(0,0,0);
                        auto f_Bj = math::make_coordinate<coordinate_type>(0,0,0);

                        const auto U_dU_attr = pot.U_dU_attr(epsilon_bp, alpha_bp, lBij, r0_bp);

                        energy += 0.5 * (1 + cos_dphi) * f1 * f2 * U_dU_attr.first;

                        // --------------------------------------------------------
                        // calc dihedral term
                        // -sin(dphi)/2 f(dtheta1) f(dtheta2) U_attr(Bij) dphi/dr
                        if(sin_dphi != real_type(0.0))
                        {
                            // remember that F = -dU. Here `coef` = -dU.
                            const auto coef = real_type(0.5) * sin_dphi *
                                              f1 * f2 * U_dU_attr.first;

                            const auto rlBij_sq = rlBij * rlBij; // 1 / |Bij|^2
                            const auto fSi = ( coef * lBij / m_lsq) * m;
                            const auto fSj = (-coef * lBij / n_lsq) * n;

                            const auto coef_Bi = dot_SBiBj * rlBij_sq;
                            const auto coef_Bj = dot_SBjBi * rlBij_sq;

                            f_Si += fSi;
                            f_Bi += (coef_Bi - real_type(1.0)) * fSi - coef_Bj * fSj;
                            f_Bj += (coef_Bj - real_type(1.0)) * fSj - coef_Bi * fSi;
                            f_Sj += fSj;
                        }

                        const auto dihd_term = real_type(0.5) * (real_type(1.0) + cos_dphi);

                        // --------------------------------------------------------
                        // calc theta1 term
                        // (1+cos(dphi))/2 df/dtheta1 f(dtheta2) U_attr(Bij) dtheta1/dr
                        if(df1 != real_type(0.0))
                        {
                            // remember that F = -dU. Here `coef` = -dU.
                            const auto coef = -dihd_term * df1 * f2 * U_dU_attr.first;

                            const auto sin_theta1 = std::sin(theta1);
                            const auto coef_rsin  = (sin_theta1 > tolerance) ?
                                       (coef / sin_theta1) : (coef / tolerance);

                            const auto fSi = (coef_rsin * rlSBi) *
                                             (cos_theta1 * BSi_reg - Bij_reg);
                            const auto fBj = (coef_rsin * rlBij) *
                                             (cos_theta1 * Bij_reg - BSi_reg);
                            f_Si += fSi;
                            f_Bi -= (fSi + fBj);
                            f_Bj += fBj;
                        }
                        // --------------------------------------------------------
                        // calc theta2 term
                        // (1+cos(dphi))/2 f(dtheta1) df/dtheta2 U_attr(Bij) dtheta2/dr
                        if(df2 != real_type(0.0))
                        {
                            // remember that F = -dU. Here `coef` = -dU.
                            const auto coef = -dihd_term * f1 * df2 * U_dU_attr.first;

                            const auto sin_theta2 = std::sin(theta2);
                            const auto coef_rsin  = (sin_theta2 > tolerance) ?
                                       (coef / sin_theta2) : (coef / tolerance);

                            const auto fBi = (coef_rsin * rlBij) *
                                             (cos_theta2 * Bji_reg - BSj_reg);
                            const auto fSj = (coef_rsin * rlSBj) *
                                             (cos_theta2 * BSj_reg - Bji_reg);
                            f_Bi += fBi;
                            f_Bj -= (fBi + fSj);
                            f_Sj += fSj;
                        }
                        // --------------------------------------------------------
                        // calc distance
                        // + 1/2 (1+cos(dphi)) f(dtheta1) f(dtheta2) dU_attr/dr  dBij/dr
                        if(U_dU_attr.second != real_type(0.0))
                        {
                            // remember that F = -dU. Here `coef` = -dU.
                            const auto coef = -dihd_term * f1 * f2 * U_dU_attr.second;
                            f_Bi += coef * Bji_reg;
                            f_Bj += coef * Bij_reg;
                        }

                        sys.force_thread(thread_id, Si) += f_Si;
                        sys.force_thread(thread_id, Bi) += f_Bi;
                        sys.force_thread(thread_id, Bj) += f_Bj;
                        sys.force_thread(thread_id, Sj) += f_Sj;

                        sys.virial_thread(thread_id)  += math::tensor_product(sys.transpose(rSi, rBi), f_Si) +
                                         math::tensor_product(                   rBi , f_Bi) +
                                         math::tensor_product(              rBi + Bij, f_Bj) +
                                         math::tensor_product(sys.transpose(rSj, rBi), f_Sj) ;
                    }
                }

                // ================================================================
                // cross stacking
                // f(theta_3) f(theta_CS) U_attr(epsilon, alpha, rij)
                //
                //       Si   Bi   Bj   Sj
                //  5'    o -- o===o -- o     3'
                //  ^    /      \ /      \    |
                //  | P o        X        o P |
                //  |    \      / \      /    v
                //  3'    o -- o===o -- o     5'
                //       Bj_next   Bj_next
                //
                // d/dr Vcs =
                //    df/dtheta3 f(theta_CS)  U_attr(eps, alp, rij) dtheta_3  /dr
                //  + f(theta_3) df/dtheta_CS U_attr(eps, alp, rij) dtheta_CS /dr
                //  + f(theta_3) f(theta_CS)  dU_attr/drij          drij/dr
                //

                const auto Bi_next        = pot.Bi_next;
                const auto Bj_next        = pot.Bj_next;
                const bool Bi_next_exists = (Bi_next != potential_type::invalid());
                const bool Bj_next_exists = (Bj_next != potential_type::invalid());

                if(!Bi_next_exists && !Bj_next_exists)
                {
                    continue; // if both interacting pair do not exist, do nothing.
                }

                // ----------------------------------------------------------------
                // calc common part, theta3 and dtheta3/dr.

                const auto cos_theta3 = math::dot_product(BSi_reg, BSj_reg);
                const auto theta3     = std::acos(math::clamp<real_type>(cos_theta3, -1, 1));
                const auto theta3_0   = parameters_.theta3_0(bp_kind);
                const auto f3         = pot.f(K_BP, pi_over_K_BP, theta3, theta3_0);
                if(f3 == real_type(0))
                {
                    // f(theta) == 0 means df(theta) is also zero.
                    // so here, both cross-stacking becomes zero. skip them.
                    continue;
                }
                const auto df3 = pot.df(K_BP, pi_over_K_BP, theta3, theta3_0);

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

                const auto K_CS            = parameters_.K_CS();
                const auto pi_over_K_CS    = parameters_.pi_over_K_CS();

                auto f_Si      = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Sj      = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Bi      = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Bj      = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Bi_next = math::make_coordinate<coordinate_type>(0,0,0);
                auto f_Bj_next = math::make_coordinate<coordinate_type>(0,0,0);

                // adjacent of Base j exists. its not the edge of the strand.
                if(Bj_next_exists)
                {
                    // ------------------------------------------------------------
                    // 5' cross stacking (the case if `Bi` is in sense strand)
                    //
                    //       Si   Bi   Bj   Sj
                    //  5'    o--> o===o <--o     3'
                    //  ^    /   `--\        \    |
                    //  | P o   tCS  \        o P |
                    //  |    \        \      /    v
                    //  3'    o -- o===o -- o     5'
                    //           Bi3   Bj5
                    //
                    // Hereafter, we use the notation illustrated above (excepting
                    // the index `Bj_next` that corresponds to Bj5).

                    const auto cs_kind = pot.cs_i_kind;

                    const auto epsilon_cs = parameters_.epsilon(cs_kind);
                    const auto alpha_cs   = parameters_.alpha(cs_kind);
                    const auto r0_cs      = parameters_.r0(cs_kind);

                    const auto& rBj5     = sys.position(Bj_next);
                    const auto  Bj5i     = sys.adjust_direction(rBj5, rBi);
                    const auto  lBj5i_sq = math::length_sq(Bj5i); // |Bj5i|^2
                    const auto  rlBj5i   = math::rsqrt(lBj5i_sq);  // 1 / |Bj5i|

                    const auto dot_theta_CS = math::dot_product(SBi, Bj5i);
                    const auto cos_theta_CS = dot_theta_CS * rlSBi * rlBj5i;
                    const auto theta_CS     =
                        std::acos(math::clamp<real_type>(cos_theta_CS, -1, 1));
                    const auto theta_CS_0   = parameters_.thetaCS_0(cs_kind);

                    const auto fCS  = pot.f(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                    // if f == 0, df is also zero. if fCS == 0, no force there
                    if(fCS != real_type(0))
                    {
                        const auto dfCS  = pot.df(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);

                        const auto lBj5i = lBj5i_sq * rlBj5i;
                        const auto U_dU_attr = pot.U_dU_attr(epsilon_cs, alpha_cs, lBj5i, r0_cs);
                        const auto Bj5i_reg  = rlBj5i * Bj5i;

                        energy += f3 * fCS * U_dU_attr.first;

                        // --------------------------------------------------------
                        // df/dtheta3 f(theta_CS)  U_attr(eps, alp, rij) dtheta_3 /dr
                        if(df3 != real_type(0))
                        {
                            const auto coef = -df3 * fCS * U_dU_attr.first;
                            f_Si += coef * fSi_theta3;
                            f_Sj += coef * fSj_theta3;
                            f_Bi += coef * fBi_theta3;
                            f_Bj += coef * fBj_theta3;
                        }
                        // --------------------------------------------------------
                        // f(theta_3) df/dtheta_CS U_attr(eps, alp, rij) dtheta_CS/dr
                        if(dfCS != real_type(0))
                        {
                            const auto coef         = -f3 * dfCS * U_dU_attr.first;
                            const auto sin_theta_CS = std::sin(theta_CS);
                            const auto coef_rsin    = (sin_theta_CS > tolerance) ?
                                       (coef / sin_theta_CS) : (coef / tolerance);

                            const auto fSi  =  coef_rsin * rlSBi  * (cos_theta_CS * BSi_reg + Bj5i_reg);
                            const auto fBj5 = -coef_rsin * rlBj5i * (cos_theta_CS * Bj5i_reg + BSi_reg);

                            f_Si      += fSi;
                            f_Bi      -= (fSi + fBj5);
                            f_Bj_next += fBj5;
                        }
                        // --------------------------------------------------------
                        // f(theta_3) f(theta_CS)  dU_attr/drij          drij/dr
                        if(U_dU_attr.second != real_type(0.0))
                        {
                            const auto coef = -f3 * fCS * U_dU_attr.second;
                            f_Bi      += coef * Bj5i_reg;
                            f_Bj_next -= coef * Bj5i_reg;
                        }
                    }
                }
                if(Bi_next_exists)
                {
                    // ------------------------------------------------------------
                    // 3' cross stacking (the case `Bi` is in a sense strand)
                    //
                    //       Si   Bi   Bj   Sj
                    //  5'    o--> o===o <--o     3'
                    //  ^    /        /--'   \    |
                    //  | P o        /  tCS   o P |
                    //  |    \      /        /    v
                    //  3'    o -- o===o -- o     5'
                    //           Bi3   Bj5
                    //
                    // Hereafter, we use the notation illustrated above (excepting
                    // the index `Bi_next` that corresponds to Bi3).

                    const auto cs_kind = pot.cs_j_kind;

                    const auto epsilon_cs = parameters_.epsilon(cs_kind);
                    const auto alpha_cs   = parameters_.alpha(cs_kind);
                    const auto r0_cs      = parameters_.r0(cs_kind);

                    const auto& rBi3     = sys.position(Bi_next);
                    const auto  Bi3j     = sys.adjust_direction(rBi3, rBj);
                    const auto  lBi3j_sq = math::length_sq(Bi3j); // |Bi3j|^2
                    const auto  rlBi3j   = math::rsqrt(lBi3j_sq);  // 1 / |Bi3j|

                    const auto dot_theta_CS = math::dot_product(SBj, Bi3j);
                    const auto cos_theta_CS = dot_theta_CS * rlSBj * rlBi3j;
                    const auto theta_CS     =
                        std::acos(math::clamp<real_type>(cos_theta_CS, -1, 1));
                    const auto theta_CS_0   = parameters_.thetaCS_0(cs_kind);

                    const auto fCS = pot.f(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                    // if f == 0, df is also zero. if fCS == 0, no force there
                    if(fCS != real_type(0))
                    {
                        const auto dfCS  = pot.df(K_CS, pi_over_K_CS, theta_CS, theta_CS_0);
                        const auto lBi3j = lBi3j_sq * rlBi3j;
                        const auto U_dU_attr = pot.U_dU_attr(epsilon_cs, alpha_cs, lBi3j, r0_cs);
                        const auto Bi3j_reg  = rlBi3j * Bi3j;

                        energy += f3 * fCS * U_dU_attr.first;

                        // --------------------------------------------------------
                        // df/dtheta3 f(theta_CS)  U_attr(eps, alp, rij) dtheta_3 /dr
                        if(df3 != real_type(0))
                        {
                            const auto coef = -df3 * fCS * U_dU_attr.first;
                            f_Si += coef * fSi_theta3;
                            f_Sj += coef * fSj_theta3;
                            f_Bi += coef * fBi_theta3;
                            f_Bj += coef * fBj_theta3;
                        }
                        // --------------------------------------------------------
                        // f(theta_3) df/dtheta_CS U_attr(eps, alp, rij) dtheta_CS/dr
                        if(dfCS != real_type(0))
                        {
                            const auto coef         = -f3 * dfCS * U_dU_attr.first;
                            const auto sin_theta_CS = std::sin(theta_CS);
                            const auto coef_rsin    = (sin_theta_CS > tolerance) ?
                                       (coef / sin_theta_CS) : (coef / tolerance);

                            const auto fSj  =  coef_rsin * rlSBj  * (cos_theta_CS * BSj_reg + Bi3j_reg);
                            const auto fBi3 = -coef_rsin * rlBi3j * (cos_theta_CS * Bi3j_reg + BSj_reg);
                            f_Sj      += fSj;
                            f_Bj      -= (fSj + fBi3);
                            f_Bi_next += fBi3;
                        }
                        // --------------------------------------------------------
                        // f(theta_3) f(theta_CS)  dU_attr/drij          drij/dr
                        if(U_dU_attr.second != real_type(0.0))
                        {
                            const auto coef = -f3 * fCS * U_dU_attr.second;
                            f_Bj      += coef * Bi3j_reg;
                            f_Bi_next -= coef * Bi3j_reg;
                        }
                    }
                }
                sys.force_thread(thread_id, Si)      += f_Si     ;
                sys.force_thread(thread_id, Sj)      += f_Sj     ;
                sys.force_thread(thread_id, Bi)      += f_Bi     ;
                sys.force_thread(thread_id, Bj)      += f_Bj     ;
                sys.force_thread(thread_id, Bi_next) += f_Bi_next;
                sys.force_thread(thread_id, Bj_next) += f_Bj_next;

                sys.virial_thread(thread_id) +=
                    math::tensor_product(sys.transpose(rSi, rBi), f_Si) +
                    math::tensor_product(                   rBi , f_Bi) +
                    math::tensor_product(              rBi + Bij, f_Bj) +
                    math::tensor_product(sys.transpose(rSj, rBi), f_Sj) +
                    math::tensor_product(sys.transpose(sys.position(Bi_next), rBi), f_Bi_next) +
                    math::tensor_product(sys.transpose(sys.position(Bj_next), rBi), f_Bj_next) ;
            }
        }
        return energy;
    }

    std::string name() const override {return "3SPN2BaseBase";}

    parameter_list_type const& parameters() const noexcept {return parameters_;}
    parameter_list_type&       parameters()       noexcept {return parameters_;}

    base_type* clone() const override
    {
        return new ThreeSPN2BaseBaseInteraction(*this);
    }

  private:

    parameter_list_type parameters_;
    partition_type partition_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize BondAngleInteraction with LocalPotentials
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif // MJOLNIR_INTERACTION_GLOBAL_3SPN_BASE_PAIRING_INTERACTION_HPP
