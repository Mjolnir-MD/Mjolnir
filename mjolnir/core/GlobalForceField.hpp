#ifndef MJOLNIR_GLOBAL_FORCE_FIELD
#define MJOLNIR_GLOBAL_FORCE_FIELD
#include "GlobalInteractionBase.hpp"
#include "GlobalPotentialBase.hpp"
#include <mjolnir/util/logger.hpp>
#include <vector>
#include <array>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class GlobalForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef GlobalPotentialBase<traitsT>      potential_base;
    typedef std::unique_ptr<potential_base>   potential_ptr;
    typedef GlobalInteractionBase<traitsT>    interaction_base;
    typedef std::unique_ptr<interaction_base> interaction_ptr;

  public:
    GlobalForceField() = default;
    ~GlobalForceField() = default;
    GlobalForceField(const GlobalForceField&) = delete;
    GlobalForceField(GlobalForceField&&) = default;
    GlobalForceField& operator=(const GlobalForceField&) = delete;
    GlobalForceField& operator=(GlobalForceField&&) = default;

    void emplace(interaction_ptr&& inter, potential_ptr&& pot)
    {
        potentials_.emplace_back(
                std::forward<interaction_ptr>(inter),
                std::forward<potential_ptr>(pot));
    }

    void      calc_force(ParticleContainer<traitsT>& pcon);
    real_type calc_energy(const ParticleContainer<traitsT>& pcon) const;

    void initialize(const ParticleContainer<traitsT>& pcon, const time_type dt);
    void reset_parameter(const std::string&, const real_type);

  private:

    struct global_potential_type
    {
        global_potential_type() = default;
        ~global_potential_type() = default;
        global_potential_type(interaction_ptr&& inter, potential_ptr&& pot)
            : interaction_(std::forward<interaction_ptr>(inter)),
              potential_(std::forward<potential_ptr>(pot))
        {}
        global_potential_type(const global_potential_type&) = delete;
        global_potential_type(global_potential_type&&) = default;
        global_potential_type& operator=(const global_potential_type&) = delete;
        global_potential_type& operator=(global_potential_type&&) = default;

        interaction_ptr interaction_;
        potential_ptr   potential_;
    };

    std::vector<global_potential_type> potentials_;

    static
    Logger& logger_;
};

template<typename traitsT>
Logger& GlobalForceField<traitsT>::logger_ =
    LoggerManager<char>::get_logger("GlobalForceField");

template<typename traitsT>
void GlobalForceField<traitsT>::initialize(
        const ParticleContainer<traitsT>& pcon, const time_type dt)
{
    for(auto iter = potentials_.begin(); iter != potentials_.end(); ++iter)
        iter->interaction_->initialize(pcon, dt);
    return;
}

template<typename traitsT>
void GlobalForceField<traitsT>::calc_force(ParticleContainer<traitsT>& pcon)
{
    for(auto iter = potentials_.begin(); iter != potentials_.end(); ++iter)
    {
        iter->interaction_->calc_force(pcon, *(iter->potential_));
        MJOLNIR_LOG_DEBUG("calculate force");
    }
    return;
}

template<typename traitsT>
typename GlobalForceField<traitsT>::real_type
GlobalForceField<traitsT>::calc_energy(
        const ParticleContainer<traitsT>& pcon) const
{
    real_type energy = 0.;
    for(auto iter = potentials_.cbegin(); iter != potentials_.cend(); ++iter)
    {
        energy += iter->interaction_->calc_energy(pcon, *(iter->potential_));
        MJOLNIR_LOG_DEBUG("calculate energy", energy);
    }
    return energy;
}

template<typename traitsT>
void GlobalForceField<traitsT>::reset_parameter(
        const std::string& name, const real_type val)
{
    for(auto iter = potentials_.begin(); iter != potentials_.end(); ++iter)
    {
        iter->potential_->reset_parameter(name, val);
    }
    return;
}


} // mjolnir

#endif /* MJOLNIR_GLOBAL_FORCE_FIELD */
