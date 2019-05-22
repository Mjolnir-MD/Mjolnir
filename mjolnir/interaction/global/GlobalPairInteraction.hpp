#ifndef MJOLNIR_INTEARACTION_GLOBAL_PAIR_INTEARACTION_HPP
#define MJOLNIR_INTEARACTION_GLOBAL_PAIR_INTEARACTION_HPP
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/string.hpp>
#include <memory>

namespace mjolnir
{

template<typename traitsT, typename potentialT, typename partitionT>
class GlobalPairInteraction final : public GlobalInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using potential_type  = potentialT;
    using partition_type  = partitionT;
    using base_type       = GlobalInteractionBase<traitsT>;
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

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

    std::string name() const override
    {return "Pair:"_s + potential_type::name();}

  private:

    potential_type potential_;
    partition_type partition_;
};

template<typename traitsT, typename potT, typename spaceT>
void GlobalPairInteraction<traitsT, potT, spaceT>::calc_force(
        system_type& sys) const noexcept
{
    for(const auto i : this->potential_.participants())
    {
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j     = ptnr.index;
            const auto& param = ptnr.parameter();

            const coordinate_type rij =
                sys.adjust_direction(sys.position(j) - sys.position(i));
            const real_type l2 = math::length_sq(rij); // |rij|^2
            const real_type rl = math::rsqrt(l2);      // 1 / |rij|
            const real_type l  = l2 * rl;              // |rij|^2 / |rij|
            const real_type f_mag = potential_.derivative(l, param);

            // if length exceeds cutoff, potential returns just 0.
            if(f_mag == 0.0){continue;}

            const coordinate_type f = rij * (f_mag * rl);
            sys.force(i) += f;
            sys.force(j) -= f;
        }
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename GlobalPairInteraction<traitsT, potT, spaceT>::real_type
GlobalPairInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type e = 0.0;
    for(const auto i : this->potential_.participants())
    {
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j     = ptnr.index;
            const auto& param = ptnr.parameter();

            const real_type l = math::length(
                sys.adjust_direction(sys.position(j) - sys.position(i)));
            e += potential_.potential(l, param);
        }
    }
    return e;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
