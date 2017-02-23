#ifndef MJOLNIR_LOCAL_INTEARACTION_BASE
#define MJOLNIR_LOCAL_INTEARACTION_BASE
#include "Particle.hpp"
#include <array>

namespace mjolnir
{

template<typename traitsT, std::size_t N>
class LocalInteractionBase;

template<typename traitsT>
class LocalInteractionBase<traitsT, 2>
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type>       particle_type;

  public:

    virtual ~LocalInteractionBase() = default;

    virtual void
    calc_force(particle_type& p1, particle_type& p2) const = 0;

    virtual real_type
    calc_energy(const particle_type& p1, const particle_type& p2) const = 0;
};

template<typename traitsT>
class LocalInteractionBase<traitsT, 3>
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type>       particle_type;

  public:

    virtual ~LocalInteractionBase() = default;

    virtual void
    calc_force(particle_type& p1, particle_type& p2, particle_type& p3) const = 0;

    virtual real_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const particle_type& p3) const = 0;
};

template<typename traitsT>
class LocalInteractionBase<traitsT, 4>
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type>       particle_type;

  public:

    virtual ~LocalInteractionBase() = default;

    virtual void
    calc_force(particle_type& p1, particle_type& p2, particle_type& p3,
               particle_type& p4) const = 0;

    virtual real_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const particle_type& p3, const particle_type& p4) const = 0;
};

} // mjolnir
#endif/* MJOLNIR_GLOBAL_INTEARACTION_BASE */
