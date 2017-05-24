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
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef System<traits_type>                   system_type;
    typedef system_type::particle_type            particle_type;
    typedef LocalInteractionBase<traitsT>         interaction_type;
    typedef std::unique_ptr<interaction_type>     interaction_ptr;
    typedef std::vector<interaction_ptr>          container_type;

  public:

    LocalForceField()  = default;
    ~LocalForceField() = default;
    LocalForceField(LocalForceField const&) = delete;
    LocalForceField(LocalForceField&&)      = default;
    LocalForceField& operator=(LocalForceField const&) = delete;
    LocalForceField& operator=(LocalForceField&&)      = default;

    void emplace(interaction_ptr&& interaction)
    {
        interactions_.emplace_back(std::move(interaction));
    }

    void      calc_force (system_type& sys);
    real_type calc_energy(const system_type& sys) const;

  private:

    container_type interactions_;
};

template<typename traitsT>
inline void LocalForceField<traitsT>::calc_force(system_type& sys)
{
    for(const auto& item : this->interactions_)
        item->calc_force(sys);
    return;
}

template<typename traitsT>
inline typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const system_type& pcon) const
{
    real_type energy = 0.0;
    for(const auto& item : this->interactions_)
        energy += item->calc_energy(sys);
    return energy;
}

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
