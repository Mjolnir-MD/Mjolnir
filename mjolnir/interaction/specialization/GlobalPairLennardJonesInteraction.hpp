#ifndef MJOLNIR_GLOBAL_PAIR_LENNARD_JONES_INTEARACTION
#define MJOLNIR_GLOBAL_PAIR_LENNARD_JONES_INTEARACTION
#include <mjolnir/interaction/GlobalPairInteraction.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <memory>

namespace mjolnir
{

// specialization for GlobalPair<LennardJones>
template<typename traitsT, typename partitionT, typename ignored_chainT>
class GlobalPairInteraction<traitsT,
    LennardJonesPotential<typename traitsT::real_type, ignored_chainT>,
    partitionT
    > final : public GlobalInteractionBase<traitsT>
{
  public:

    using traits_type     = traitsT;
    using partition_type  = partitionT;
    using base_type       = GlobalInteractionBase<traits_type>;
    using potential_type  =
        LennardJonesPotential<typename traitsT::real_type, ignored_chainT>;

    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

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
        MJOLNIR_SCOPE(GlobalPairInteraction<LennardJonesPotential>::initialize(), 0);
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->partition_.initialize(sys, this->potential_);
        this->partition_.update(sys, this->potential_);
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(GlobalPairInteraction<LennardJonesPotential>::update(), 0);
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.update(sys);
        // potential update may change the cutoff length!
        this->partition_.reconstruct(sys, this->potential_);
    }

    void calc_force(system_type& sys) override
    {
        partition_.update(sys, this->potential_);

        constexpr auto  cutoff_ratio    = potential_type::cutoff_ratio;
        constexpr auto  cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;
        const     auto& param           = potential_.parameters();
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j     = ptnr.index;
                const auto& param = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys[j].position - sys[i].position);
                const real_type l_sq = length_sq(rij);

                const real_type sigma_sq = param.first * param.first;
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                const real_type epsilon = param.second;

                const real_type rcp_l_sq = 1 / l_sq;
                const real_type s2l2 = sigma_sq * rcp_l_sq;
                const real_type s6l6 = s2l2 * s2l2 * s2l2;

                const coordinate_type f = rij *
                    (24 * epsilon * (s6l6 - 2 * s6l6 * s6l6) * rcp_l_sq);

                sys[i].force += f;
                sys[j].force -= f;
            }
        }
        return ;
    }

    real_type calc_energy(const system_type& sys) const override
    {
        real_type E(0);

        constexpr auto  cutoff_ratio    = potential_type::cutoff_ratio;
        constexpr auto  cutoff_ratio_sq = cutoff_ratio * cutoff_ratio;
        constexpr auto  coef_at_cutoff  = potential_type::coef_at_cutoff;
        const     auto& param           = potential_.parameters();
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j     = ptnr.index;
                const auto& param = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys[j].position - sys[i].position);
                const real_type l_sq = length_sq(rij);

                const real_type sigma_sq = param.first * param.first;
                if(sigma_sq * cutoff_ratio_sq < l_sq) {continue;}

                const real_type epsilon = param.second;

                const real_type s2l2 = sigma_sq / l_sq;
                const real_type s6l6 = s2l2 * s2l2 * s2l2;

                E += 4 * epsilon * (s6l6 * s6l6 - s6l6 - coef_at_cutoff);
            }
        }
        return E;
    }

    std::string name() const override
    {return "GlobalPairLennardJonesInteraction";}

  private:

    potential_type potential_;
    partition_type partition_;
};

} // mjolnir
#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
