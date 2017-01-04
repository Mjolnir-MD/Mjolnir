#ifndef MJOLNIR_GLOBAL_INTEARACTION_BASE
#define MJOLNIR_GLOBAL_INTEARACTION_BASE
#include "ParticleContainer.hpp" 
#include "GlobalPotentialBase.hpp"

namespace mjolnir
{

template<typename traitsT>
class GlobalInteractionBase
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef GlobalPotentialBase<traits_type> potential_type;

  public:
    virtual ~GlobalInteractionBase() = default;

    virtual void
    calc_force(particle_container_type& pcon, potential_type& pot) = 0;

    virtual real_type
    calc_energy(const particle_container_type& pcon,
                const potential_type& pot) const = 0;

    virtual void initialize(
            const particle_container_type& pcon, const time_type dt) = 0;
};

} // mjolnir
#endif/* MJOLNIR_GLOBAL_INTEARACTION_BASE */
