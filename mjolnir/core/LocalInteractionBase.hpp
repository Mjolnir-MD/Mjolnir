#ifndef MJOLNIR_LOCAL_INTEARACTION_BASE
#define MJOLNIR_LOCAL_INTEARACTION_BASE
#include "Particle.hpp"
#include "LocalPotentialBase.hpp"
#include <array>

namespace mjolnir
{

template<typename traitsT, std::size_t N>
class LocalInteractionBase
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type>       particle_type;
    typedef particle_type*                  particle_ptr; // XXX
    typedef std::array<particle_ptr, N>     particle_ptrs;
    typedef LocalPotentialBase<traits_type> potential_type;

  public:

    virtual ~LocalInteractionBase() = default;

    virtual void
    calc_force(particle_ptrs ps, const potential_type& pot) const = 0;

    virtual real_type
    calc_energy(const particle_ptrs ps, const potential_type& pot) const = 0;
};

} // mjolnir
#endif/* MJOLNIR_GLOBAL_INTEARACTION_BASE */
