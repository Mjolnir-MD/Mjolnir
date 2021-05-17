#ifndef MJOLNIR_FORCEFIELD_3SPN2_BASE_STACKING_INTEARACTION_HPP
#define MJOLNIR_FORCEFIELD_3SPN2_BASE_STACKING_INTEARACTION_HPP
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseStackingPotential.hpp>
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/string.hpp>
#include <memory>

namespace mjolnir
{

// It calculates a stacking energy/force that is a part of 3SPN2 DNA model.
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
// And also 3SPN2.C model.
// - G. S. Freeman, D. M. Hinckley, J. P. Lequieu, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2014)
//
// The potential function closely tied to the interaction, so this interaction
// class is implemented with its own, uninterchangeable potential class.
// It does not take any potential class as a template parameter because
// exchanging potential function of this interaction does not make any sense.
//
// Note: an identifier starts with a digit is not allowed in C++ standard.
//       see N3337 2.11 for detail. So `3SPN2BaseStacking` is not a valid name.
template<typename traitsT>
class ThreeSPN2BaseStackingInteraction final : public LocalInteractionBase<traitsT>
{
  public:
    using traits_type          = traitsT;
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;

    using potential_type       = ThreeSPN2BaseStackingPotential<real_type>;
    using base_kind            = parameter_3SPN2::base_kind;
    using base_stack_kind      = parameter_3SPN2::base_stack_kind;

    using indices_type         = std::array<std::size_t, 3>;
    using parameter_type       = base_stack_kind;
    using parameter_index_pair = std::pair<indices_type, parameter_type>;
    using container_type       = std::vector<parameter_index_pair>;

    // to register adjacent nucleotides to Topology...
    using nucleotide_index_type = parameter_3SPN2::NucleotideIndex;

  public:

    ThreeSPN2BaseStackingInteraction(const connection_kind_type kind,
            container_type&& para, potential_type&& pot,
            std::vector<nucleotide_index_type>&& nuc_idx)
        : kind_(kind), parameters_(std::move(para)), potential_(std::move(pot)),
          nucleotide_index_(std::move(nuc_idx))
    {}
    ~ThreeSPN2BaseStackingInteraction() {}
    ThreeSPN2BaseStackingInteraction(const ThreeSPN2BaseStackingInteraction&) = default;
    ThreeSPN2BaseStackingInteraction(ThreeSPN2BaseStackingInteraction&&)      = default;
    ThreeSPN2BaseStackingInteraction& operator=(const ThreeSPN2BaseStackingInteraction&) = default;
    ThreeSPN2BaseStackingInteraction& operator=(ThreeSPN2BaseStackingInteraction&&)      = default;

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        this->potential_.initialize(sys);
    }

    void update(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        this->potential_.update(sys);
    }

    // do nothing. this is used to reduce margin of neighbor list, and added
    // to this class for the consistency.
    void reduce_margin(const real_type, const system_type&) override {return;}
    void  scale_margin(const real_type, const system_type&) override {return;}

    void      calc_force (system_type&)           const noexcept override;
    real_type calc_energy(const system_type&)     const noexcept override;
    real_type calc_force_and_energy(system_type&) const noexcept override;

    std::string name() const override {return "3SPN2BaseStacking"_s;}

    // Unlike other interactions, it registers edges between adjacent nucleotides.
    // Because Base-Base interaction counts number of nucleotides that separates
    // particles.
    void write_topology(topology_type& topol) const override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        if(this->kind_.empty() || this->kind_ == "none")
        {
            MJOLNIR_LOG_WARN("3SPN2 Base-Base Interaction (base pairing + "
                             "cross stacking) requires the number of nucleotides"
                             " that separates bases but topology is not set.");
            MJOLNIR_LOG_WARN("I trust that you know what you are doing.");
            return;
        }

        for(std::size_t i=1; i<nucleotide_index_.size(); ++i)
        {
            constexpr auto nil = nucleotide_index_type::nil();
            const auto& Base5 = nucleotide_index_.at(i-1);
            const auto& Base3 = nucleotide_index_.at(i);

            if(Base5.strand != Base3.strand) {continue;}

            topol.add_connection(Base5.S, Base3.S, this->kind_);
            topol.add_connection(Base5.S, Base3.B, this->kind_);
            topol.add_connection(Base5.B, Base3.S, this->kind_);
            topol.add_connection(Base5.B, Base3.B, this->kind_);

            if(Base5.P != nil)
            {
                topol.add_connection(Base5.P, Base3.S, this->kind_);
                topol.add_connection(Base5.P, Base3.B, this->kind_);
            }
            if(Base3.P != nil)
            {
                topol.add_connection(Base5.S, Base3.P, this->kind_);
                topol.add_connection(Base5.B, Base3.P, this->kind_);
            }
            if(Base5.P != nil && Base3.P != nil)
            {
                topol.add_connection(Base5.P, Base3.P, this->kind_);
            }
        }
        return;
    }

    connection_kind_type const& connection_kind() const noexcept {return kind_;}

    container_type const& parameters() const noexcept {return parameters_;}
    container_type&       parameters()       noexcept {return parameters_;}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type&       potential()       noexcept {return potential_;}

    base_type* clone() const override
    {
        return new ThreeSPN2BaseStackingInteraction(kind_,
                container_type(parameters_), potential_type(potential_),
                std::vector<nucleotide_index_type>(nucleotide_index_));
    }

  private:

    connection_kind_type kind_;
    container_type parameters_;
    potential_type potential_;
    std::vector<nucleotide_index_type> nucleotide_index_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT>
void ThreeSPN2BaseStackingInteraction<traitsT>::calc_force(
        system_type& sys) const noexcept
{
    constexpr auto tolerance = math::abs_tolerance<real_type>();
    for(const auto& idxp : this->parameters_)
    {
        // ====================================================================
        // Base Stacking
        // U_BS = U_rep(rij) + f(theta) U_attr(rij)
        //
        //        SBi
        //     Si --> Bi
        //    /     `-^
        //   Pj theta | rij
        //    \       |
        //     Sj --- Bj
        //
        // dU_BS/dr = dU_rep(rij)/dr + df(theta) U_attr(rij) dtheta/dr
        //                           + f(theta) dU_attr(rij) drij/dr

        const std::size_t Si = idxp.first[0];
        const std::size_t Bi = idxp.first[1];
        const std::size_t Bj = idxp.first[2];
        const auto   bs_kind = idxp.second;

        const auto& rBi = sys.position(Bi);
        const auto& rBj = sys.position(Bj);
        const auto& rSi = sys.position(Si);

        const auto Bji = sys.adjust_direction(rBj, rBi); // Bj -> Bi
        const auto SBi = sys.adjust_direction(rSi, rBi); // Si -> Bi

        const auto lBji_sq = math::length_sq(Bji); // |Bji|^2
        const auto rlBji   = math::rsqrt(lBji_sq); // 1 / |Bji|
        const auto lBji    = lBji_sq * rlBji;      // |Bji|
        const auto Bji_reg = Bji * rlBji;          // Bji / |Bji|

        auto f_Bi = math::make_coordinate<coordinate_type>(0,0,0);
        auto f_Bj = math::make_coordinate<coordinate_type>(0,0,0);
        auto f_Si = math::make_coordinate<coordinate_type>(0,0,0);

        // ====================================================================
        // calc repulsive part, which does not depend on angle term.

        const auto dU_rep = potential_.dU_rep(bs_kind, lBji);
        if(dU_rep != real_type(0.0))
        {
            f_Bi -= dU_rep * Bji_reg;
            f_Bj += dU_rep * Bji_reg;
        }

        // --------------------------------------------------------------------
        // calc theta
        //
        //        SBi
        //     Si --> Bi
        //    /     `-^
        //   Pj theta | rij
        //    \       |
        //     Sj --- Bj

        const auto lSBi_sq = math::length_sq(SBi);
        const auto rlSBi   = math::rsqrt(lSBi_sq);
        const auto SBi_reg = SBi * rlSBi;

        const auto cos_theta = math::dot_product(SBi_reg, Bji_reg);
        const auto theta = std::acos(math::clamp<real_type>(cos_theta, -1, 1));

        // ====================================================================
        // calc attractive part
        //
        // dU_BS^attr/dr = df(theta) U_attr(rij) dtheta/dr
        //               + f(theta) dU_attr(rij) drij/dr

        const auto theta_0 = potential_.theta_0(bs_kind);
        const auto f_theta = potential_.f(theta, theta_0);
        if(f_theta == real_type(0.0))
        {
            // purely repulsive. add the repulsive force to the system and quit.
            sys.force(Bi) += f_Bi;
            sys.force(Bj) += f_Bj;
            // Bji = Bj -> Bi = Bi - Bj
            sys.virial() += math::tensor_product(Bji, f_Bi); // (Bi - Bj) * fBi
            continue;
        }
        const auto df_theta  = potential_.df(theta, theta_0);
        const auto U_dU_attr = potential_.U_dU_attr(bs_kind, lBji);

        // --------------------------------------------------------------------
        // calc the first term in the attractive part
        // = df(theta) U_attr(rij) dtheta/dr
        //
        if(df_theta != real_type(0.0))
        {
            const auto coef = -df_theta * U_dU_attr.first;

            const auto sin_theta = std::sin(theta);
            const auto coef_rsin = (sin_theta > tolerance) ?
                                   coef / sin_theta : coef / tolerance;

            const auto fSi = (coef_rsin * rlSBi) * (Bji_reg - cos_theta * SBi_reg);
            const auto fBj = (coef_rsin * rlBji) * (SBi_reg - cos_theta * Bji_reg);

            f_Si +=  fSi;
            f_Bi -= (fSi + fBj);
            f_Bj +=        fBj;
        }

        // --------------------------------------------------------------------
        // calc the second term in the attractive part
        // = f(theta) dU_attr(rij) drij/dr

        if(U_dU_attr.second != real_type(0.0))
        {
            const auto coef = f_theta * U_dU_attr.second;
            f_Bi -= coef * Bji_reg;
            f_Bj += coef * Bji_reg;
        }

        sys.force(Bi) += f_Bi;
        sys.force(Bj) += f_Bj;
        sys.force(Si) += f_Si;

        sys.virial() += math::tensor_product(rBi,       f_Bi) +
                        math::tensor_product(rBi - Bji, f_Bj) +
                        math::tensor_product(rBi - SBi, f_Si);
    }
    return ;
}

template<typename traitsT>
typename ThreeSPN2BaseStackingInteraction<traitsT>::real_type
ThreeSPN2BaseStackingInteraction<traitsT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const auto& idxp : this->parameters_)
    {
        // ====================================================================
        // Base Stacking
        //
        //        SBi
        //     Si --> Bi
        //    /     `-^
        //   Pj    th | rij
        //    \       |
        //     Sj --- Bj
        //
        // U_BS = U_rep(rij) + f(theta) U_attr(rij)

        const std::size_t Si = idxp.first[0];
        const std::size_t Bi = idxp.first[1];
        const std::size_t Bj = idxp.first[2];
        const auto   bs_kind = idxp.second;

        const auto& rBi = sys.position(Bi);
        const auto& rBj = sys.position(Bj);
        const auto& rSi = sys.position(Si);

        const auto Bji = sys.adjust_direction(rBj, rBi); // Bj -> Bi
        const auto SBi = sys.adjust_direction(rSi, rBi); // Si -> Bi

        const auto lBji_sq = math::length_sq(Bji);
        const auto rlBji   = math::rsqrt(lBji_sq);
        const auto lBji    = lBji_sq * rlBji;

        // ====================================================================
        // calc repulsive part, which does not depend on angle term.

        E += potential_.U_rep(bs_kind, lBji);

        // --------------------------------------------------------------------
        // calc theta

        const auto lSBi_sq = math::length_sq(SBi);
        const auto rlSBi   = math::rsqrt(lSBi_sq);

        const auto dot_SBij = math::dot_product(SBi, Bji);
        const auto cos_SBij = dot_SBij * rlBji * rlSBi;
        const auto theta = std::acos(math::clamp<real_type>(cos_SBij, -1, 1));

        // ====================================================================
        // calc attractive part
        // = f(theta) U_attr(rij)

        const auto theta_0 = potential_.theta_0(bs_kind);
        const auto f_theta = potential_.f(theta, theta_0);
        if(f_theta == real_type(0.0))
        {
            continue; // purely repulsive.
        }
        const auto U_attr = potential_.U_attr(bs_kind, lBji);
        E += U_attr * f_theta;
    }
    return E;
}

template<typename traitsT>
typename ThreeSPN2BaseStackingInteraction<traitsT>::real_type
ThreeSPN2BaseStackingInteraction<traitsT>::calc_force_and_energy(
        system_type& sys) const noexcept
{
    constexpr auto tolerance = math::abs_tolerance<real_type>();
    real_type energy = 0;
    for(const auto& idxp : this->parameters_)
    {
        // ====================================================================
        // Base Stacking
        // U_BS = U_rep(rij) + f(theta) U_attr(rij)
        //
        //        SBi
        //     Si --> Bi
        //    /     `-^
        //   Pj theta | rij
        //    \       |
        //     Sj --- Bj
        //
        // dU_BS/dr = dU_rep(rij)/dr + df(theta) U_attr(rij) dtheta/dr
        //                           + f(theta) dU_attr(rij) drij/dr

        const std::size_t Si = idxp.first[0];
        const std::size_t Bi = idxp.first[1];
        const std::size_t Bj = idxp.first[2];
        const auto   bs_kind = idxp.second;

        const auto& rBi = sys.position(Bi);
        const auto& rBj = sys.position(Bj);
        const auto& rSi = sys.position(Si);

        const auto Bji = sys.adjust_direction(rBj, rBi); // Bj -> Bi
        const auto SBi = sys.adjust_direction(rSi, rBi); // Si -> Bi

        const auto lBji_sq = math::length_sq(Bji); // |Bji|^2
        const auto rlBji   = math::rsqrt(lBji_sq); // 1 / |Bji|
        const auto lBji    = lBji_sq * rlBji;      // |Bji|
        const auto Bji_reg = Bji * rlBji;          // Bji / |Bji|

        auto f_Bi = math::make_coordinate<coordinate_type>(0,0,0);
        auto f_Bj = math::make_coordinate<coordinate_type>(0,0,0);
        auto f_Si = math::make_coordinate<coordinate_type>(0,0,0);

        // ====================================================================
        // calc repulsive part, which does not depend on angle term.

        const auto dU_rep = potential_.dU_rep(bs_kind, lBji);
        if(dU_rep != real_type(0.0))
        {
            f_Bi -= dU_rep * Bji_reg;
            f_Bj += dU_rep * Bji_reg;
        }
        energy += potential_.U_rep(bs_kind, lBji);

        // --------------------------------------------------------------------
        // calc theta
        //
        //        SBi
        //     Si --> Bi
        //    /     `-^
        //   Pj theta | rij
        //    \       |
        //     Sj --- Bj

        const auto lSBi_sq = math::length_sq(SBi);
        const auto rlSBi   = math::rsqrt(lSBi_sq);
        const auto SBi_reg = SBi * rlSBi;

        const auto cos_theta = math::dot_product(SBi_reg, Bji_reg);
        const auto theta = std::acos(math::clamp<real_type>(cos_theta, -1, 1));

        // ====================================================================
        // calc attractive part
        //
        // dU_BS^attr/dr = df(theta) U_attr(rij) dtheta/dr
        //               + f(theta) dU_attr(rij) drij/dr

        const auto theta_0 = potential_.theta_0(bs_kind);
        const auto f_theta = potential_.f(theta, theta_0);
        if(f_theta == real_type(0.0))
        {
            // purely repulsive. add the repulsive force to the system and quit.
            sys.force(Bi) += f_Bi;
            sys.force(Bj) += f_Bj;
            sys.virial() += math::tensor_product(Bji, f_Bi); // (Bi - Bj) * f_Bi
            continue;
        }
        const auto df_theta  = potential_.df(theta, theta_0);
        const auto U_dU_attr = potential_.U_dU_attr(bs_kind, lBji);

        energy += U_dU_attr.first * f_theta;

        // --------------------------------------------------------------------
        // calc the first term in the attractive part
        // = df(theta) U_attr(rij) dtheta/dr
        //
        if(df_theta != real_type(0.0))
        {
            const auto coef = -df_theta * U_dU_attr.first;

            const auto sin_theta = std::sin(theta);
            const auto coef_rsin = (sin_theta > tolerance) ?
                                   coef / sin_theta : coef / tolerance;

            const auto fSi = (coef_rsin * rlSBi) * (Bji_reg - cos_theta * SBi_reg);
            const auto fBj = (coef_rsin * rlBji) * (SBi_reg - cos_theta * Bji_reg);

            f_Si +=  fSi;
            f_Bi -= (fSi + fBj);
            f_Bj +=        fBj;
        }

        // --------------------------------------------------------------------
        // calc the second term in the attractive part
        // = f(theta) dU_attr(rij) drij/dr

        if(U_dU_attr.second != real_type(0.0))
        {
            const auto coef = f_theta * U_dU_attr.second;
            f_Bi -= coef * Bji_reg;
            f_Bj += coef * Bji_reg;
        }

        sys.force(Bi) += f_Bi;
        sys.force(Bj) += f_Bj;
        sys.force(Si) += f_Si;

        sys.virial() += math::tensor_product(rBi,       f_Bi) +
                        math::tensor_product(rBi - Bji, f_Bj) +
                        math::tensor_product(rBi - SBi, f_Si);
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

extern template class ThreeSPN2BaseStackingInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ThreeSPN2BaseStackingInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ThreeSPN2BaseStackingInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ThreeSPN2BaseStackingInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD
#endif // MJOLNIR_INTEARACTION_GLOBAL_3SPN_BASE_STACKING_INTEARACTION_HPP
