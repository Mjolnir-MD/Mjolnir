#ifndef MJOLNIR_LOCAL_FORCE_FIELD
#define MJOLNIR_LOCAL_FORCE_FIELD
#include "LocalInteractionBase.hpp"
#include <mjolnir/util/logger.hpp>
#include <utility>
#include <vector>
#include <array>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class LocalForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type       time_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type>             particle_type;
    typedef ParticleContainer<traits_type>        particle_container_type;
    typedef LocalInteractionBase<traitsT>         interaction_type;
    typedef std::unique_ptr<interaction_type>     interaction_ptr;
    typedef std::vector<interaction_ptr>          container_type;

  public:

    LocalForceField() = default;
    ~LocalForceField() = default;
    LocalForceField(LocalForceField const&) = delete;
    LocalForceField(LocalForceField&&)      = default;
    LocalForceField& operator=(LocalForceField const&) = delete;
    LocalForceField& operator=(LocalForceField&&)      = default;

    void emplace(interaction_ptr&& interaction);

    void      calc_force(particle_container_type& pcon);
    real_type calc_energy(const particle_container_type& pcon) const;

    void reset_parameter(const std::string&, const real_type);

  private:

    container_type interactions_;
    static
    Logger& logger_;
};

template<typename traitsT>
Logger& LocalForceField<traitsT>::logger_ =
    LoggerManager<char>::get_logger("LocalForceField");

template<typename traitsT>
inline void LocalForceField<traitsT>::emplace(
        interaction_ptr&& interaction)
{
    interactions_.emplace_back(std::forward<interaction_ptr>(interaction));
    return;
}

template<typename traitsT>
void LocalForceField<traitsT>::calc_force(particle_container_type& pcon)
{
    for(auto iter = interactions_.cbegin(); iter != interactions_.cend(); ++iter)
        (*iter)->calc_force(pcon);
    return;
}

template<typename traitsT>
typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const particle_container_type& pcon) const
{
    real_type energy = 0.0;
    for(auto iter = interactions_.cbegin(); iter != interactions_.cend(); ++iter)
    {
        const real_type e = (*iter)->calc_energy(pcon);
        MJOLNIR_LOG_DEBUG("calculate energy", e);
        energy += e;
    }
    return energy;
}

template<typename traitsT>
void LocalForceField<traitsT>::reset_parameter(
        const std::string& name, const real_type val)
{
    for(auto iter = interactions_.begin(); iter != interactions_.end(); ++iter)
        (*iter)->reset_parameter(name, val);
    return;
}

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
