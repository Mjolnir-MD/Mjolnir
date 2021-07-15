#ifndef MJOLNIR_INTERACTION_GLOBAL_PAIR_EXCLUDED_VOLUME_INTEARACTION_HPP
#define MJOLNIR_INTERACTION_GLOBAL_PAIR_EXCLUDED_VOLUME_INTEARACTION_HPP
#include <mjolnir/forcefield/global/GlobalPairInteraction.hpp>
#include <mjolnir/forcefield/global/ExcludedVolumePotential.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <memory>

namespace mjolnir
{

// It is a specialization of GlobalPairInteraction for ExcludedVolumePotential.
// In the case of ExcludedVolumePotential, we can omit `sqrt` call that is
// normally used to calculate distance because we only needs the squared distance.
template<typename realT, template<typename, typename> class boundaryT>
class GlobalPairInteraction<
    SimulatorTraits<realT, boundaryT>,
    ExcludedVolumePotential<realT>
    > final : public GlobalInteractionBase<SimulatorTraits<realT, boundaryT>>
{
  public:

    using traits_type     = SimulatorTraits<realT, boundaryT>;
    using base_type       = GlobalInteractionBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using topology_type   = typename base_type::topology_type;
    using boundary_type   = typename base_type::boundary_type;
    using potential_type  = ExcludedVolumePotential<real_type>;
    using partition_type  = SpatialPartition<traits_type, potential_type>;
    using parameter_list_type = ParameterList<traits_type, potential_type>;

  public:

    GlobalPairInteraction(potential_type&& pot, parameter_list_type&& para,
                          partition_type&& part)
        : potential_(std::move(pot)), parameters_(std::move(para)),
          partition_(std::move(part))
    {}
    ~GlobalPairInteraction() override {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());

        this->potential_ .initialize(sys);
        this->parameters_.initialize(sys, topol, potential_);
        this->partition_ .initialize(sys, parameters_.cref());
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());

        this->potential_ .update(sys);
        this->parameters_.update(sys, topol, potential_);
        this->partition_.initialize(sys, this->parameters_.cref());
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->partition_.reduce_margin(dmargin, sys, this->parameters_.cref());
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->partition_.scale_margin(scale, sys, this->parameters_.cref());
        return;
    }

    void calc_force(system_type& sys) const noexcept override
    {
        const auto cutoff_ratio    = potential_.cutoff_ratio();
        const auto cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;

        const auto epsilon12 = potential_.epsilon() * 12;

        const auto leading_participants = this->parameters_.leading_participants();
        for(const std::size_t i : leading_participants)
        {
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j    = ptnr.index;
                const auto& para = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l_sq = math::length_sq(rij);

                const real_type sigma_sq = para.radius * para.radius;
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                const real_type rcp_l_sq = real_type(1) / l_sq;
                const real_type s2l2     = sigma_sq * rcp_l_sq;
                const real_type s6l6     = s2l2 * s2l2 * s2l2;

                const coordinate_type f = rij *
                    (-epsilon12 * s6l6 * s6l6 * rcp_l_sq);

                sys.force(i) += f;
                sys.force(j) -= f;

                // rij * Fj = (rj - ri) * Fj = (ri - rj) * Fi
                sys.virial() += math::tensor_product(rij, -f);
            }
        }
        return ;
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type energy(0);

        const auto coef_at_cutoff  = potential_.coef_at_cutoff();
        const auto cutoff_ratio    = potential_.cutoff_ratio();
        const auto cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;

        const auto epsilon = potential_.epsilon();

        const auto leading_participants = this->parameters_.leading_participants();
        for(const std::size_t i : leading_participants)
        {
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j    = ptnr.index;
                const auto& para = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l_sq = math::length_sq(rij);

                const real_type sigma_sq = para.radius * para.radius;
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                const real_type s2l2 = sigma_sq / l_sq;
                const real_type s6l6 = s2l2 * s2l2 * s2l2;

                energy += epsilon * (s6l6 * s6l6 - coef_at_cutoff);
            }
        }
        return energy;
    }

    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        const auto coef_at_cutoff  = potential_.coef_at_cutoff();
        const auto cutoff_ratio    = potential_.cutoff_ratio();
        const auto cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;

        const auto epsilon   = potential_.epsilon();
        const auto epsilon12 = epsilon * 12;

        real_type energy = 0;
        const auto leading_participants = this->parameters_.leading_participants();
        for(const std::size_t i : leading_participants)
        {
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j    = ptnr.index;
                const auto& para = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l_sq = math::length_sq(rij);

                const real_type sigma_sq = para.radius * para.radius;
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                const real_type rcp_l_sq = real_type(1) / l_sq;
                const real_type s2l2     = sigma_sq * rcp_l_sq;
                const real_type s6l6     = s2l2 * s2l2 * s2l2;

                energy += epsilon * (s6l6 * s6l6 - coef_at_cutoff);

                const auto f = rij * (-epsilon12 * s6l6 * s6l6 * rcp_l_sq);

                sys.force(i) += f;
                sys.force(j) -= f;

                // rij * Fj = (rj - ri) * Fj = (ri - rj) * Fi
                sys.virial() += math::tensor_product(rij, -f);
            }
        }
        return energy;
    }

    std::string name() const override {return "GlobalPairExcludedVolume";}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type&       potential()       noexcept {return potential_;}

    partition_type const& partition() const noexcept {return partition_;}

    base_type* clone() const override
    {
        return new GlobalPairInteraction(*this);
    }

  private:

    potential_type      potential_;
    parameter_list_type parameters_;
    partition_type      partition_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{
// ============================================================================
// exv
// CellList
extern template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<double>>;
extern template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<float >>;
extern template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>>;
extern template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
