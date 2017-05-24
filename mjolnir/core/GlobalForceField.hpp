#ifndef MJOLNIR_GLOBAL_FORCE_FIELD
#define MJOLNIR_GLOBAL_FORCE_FIELD
#include "GlobalInteractionBase.hpp"
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
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename system_type::particle_type   particle_type;
    typedef GlobalInteractionBase<traitsT>    interaction_base;
    typedef std::unique_ptr<interaction_base> interaction_ptr;

  public:
    GlobalForceField() = default;
    ~GlobalForceField() = default;
    GlobalForceField(const GlobalForceField&) = delete;
    GlobalForceField(GlobalForceField&&)      = default;
    GlobalForceField& operator=(const GlobalForceField&) = delete;
    GlobalForceField& operator=(GlobalForceField&&)      = default;

    void emplace(interaction_ptr&& inter)
    {
        interactions_.emplace_back(std::move(inter));
    }

    void initialize(const system_type& sys, const real_type dt)
    {
        for(auto& item : this->interactions_)
            item->initialize(sys, dt);
    }

    void      calc_force (system_type&)       const;
    real_type calc_energy(const system_type&) const;

  private:

    std::vector<interaction_ptr> interactions_;

    static
    Logger& logger_;
};

template<typename traitsT>
Logger& GlobalForceField<traitsT>::logger_ =
    LoggerManager<char>::get_logger("GlobalForceField");

template<typename traitsT>
inline void GlobalForceField<traitsT>::calc_force(system_type& sys) const
{
    for(const auto& item : this->interactions_)
        item->calc_force(sys);
    return;
}

template<typename traitsT>
inline typename GlobalForceField<traitsT>::real_type
GlobalForceField<traitsT>::calc_energy(const system_type& pcon) const
{
    real_type energy = 0.;
    for(const auto& item : this->interactions_)
        energy += item->calc_energy(pcon);
    return energy;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_FORCE_FIELD */
