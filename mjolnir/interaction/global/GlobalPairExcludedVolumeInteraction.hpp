#ifndef MJOLNIR_INTERACTION_GLOBAL_PAIR_EXCLUDED_VOLUME_INTEARACTION_HPP
#define MJOLNIR_INTERACTION_GLOBAL_PAIR_EXCLUDED_VOLUME_INTEARACTION_HPP
#include <mjolnir/interaction/global/GlobalPairInteraction.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
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
    using boundary_type   = typename base_type::boundary_type;
    using potential_type  = ExcludedVolumePotential<real_type>;
    using partition_type  = SpatialPartition<traits_type, potential_type>;

  public:
    GlobalPairInteraction()  = default;
    ~GlobalPairInteraction() = default;

    GlobalPairInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.initialize(sys);
        this->partition_.initialize(sys, this->potential_);
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.update(sys);
        // potential update may change the cutoff length!
        this->partition_.initialize(sys, this->potential_);
    }

    void update_margin(const real_type dmargin, const system_type& sys) override
    {
        this->partition_.update(dmargin, sys, this->potential_);
        return;
    }

    void calc_force(system_type& sys) const noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();

        const auto cutoff_ratio    = potential_.cutoff_ratio();
        const auto cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;
        const auto epsilon12       = 12 * potential_.epsilon();
        for(const auto i : this->potential_.participants())
        {
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j     = ptnr.index;
                const auto& param = ptnr.parameter(); // sum of radius

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(j) - sys.position(i));
                const real_type l_sq = math::length_sq(rij);

                const real_type sigma_sq = param * param;
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                MJOLNIR_LOG_DEBUG("calculating force between ", i, " and ", j);

                const real_type rcp_l_sq = real_type(1) / l_sq;
                const real_type s2l2     = sigma_sq * rcp_l_sq;
                const real_type s6l6     = s2l2 * s2l2 * s2l2;

                const coordinate_type f = rij *
                    (-epsilon12 * s6l6 * s6l6 * rcp_l_sq);

                sys.force(i) += f;
                sys.force(j) -= f;
            }
        }
        return ;
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();

        real_type E(0);

        const auto cutoff_ratio    = potential_.cutoff_ratio();
        const auto cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;
        const auto coef_at_cutoff  = potential_.coef_at_cutoff();
        const auto epsilon         = potential_.epsilon();
        for(const auto i : this->potential_.participants())
        {
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j     = ptnr.index;
                const auto& param = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(j) - sys.position(i));
                const real_type l_sq = math::length_sq(rij);

                const real_type sigma_sq = param * param;
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                MJOLNIR_LOG_DEBUG("calculating energy between ", i, " and ", j);

                const real_type s2l2 = sigma_sq / l_sq;
                const real_type s6l6 = s2l2 * s2l2 * s2l2;

                E += epsilon * (s6l6 * s6l6 - coef_at_cutoff);
            }
        }
        return E;
    }

    std::string name() const override
    {return "GlobalPairExcludedVolumeInteraction";}

  private:

    potential_type potential_;
    partition_type partition_;
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
extern template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>>;
extern template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float> >;
extern template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>>;
extern template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float> >;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
