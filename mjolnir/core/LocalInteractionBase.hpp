#ifndef MJOLNIR_LOCAL_INTEARACTION_BASE
#define MJOLNIR_LOCAL_INTEARACTION_BASE
#include "ParticleContainer.hpp"
#include <array>

namespace mjolnir
{

template<typename traitsT>
class LocalInteractionBase
{
  public:

    typedef traitsT traits_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename particle_container_type::particle_type particle_type;

  public:

    virtual ~LocalInteractionBase() = default;

    virtual void
    calc_force(particle_container_type& pcon) const = 0;

    virtual real_type
    calc_energy(const particle_container_type& pcon) const = 0;

    virtual void
    reset_parameter(const std::string&, const real_type) = 0;
};

} // mjolnir
#endif/* MJOLNIR_GLOBAL_INTEARACTION_BASE */
