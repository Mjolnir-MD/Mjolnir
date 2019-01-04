#ifndef MJOLNIR_EXTERNAL_FORCE_INTEARACTION_BASE
#define MJOLNIR_EXTERNAL_FORCE_INTEARACTION_BASE
#include <mjolnir/core/System.hpp>
#include <string>

namespace mjolnir
{

template<typename traitsT>
class ExternalForceInteractionBase
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename traits_type::boundary_type   boundary_type;

  public:

    virtual ~ExternalForceInteractionBase() = default;

    virtual void initialize(const system_type&) = 0;
    virtual void update    (const system_type&) = 0;

    virtual void      calc_force (system_type&)             = 0;
    virtual real_type calc_energy(const system_type&) const = 0;

    virtual std::string name() const = 0;
};

} // mjolnir
#endif/* MJOLNIR_EXTERNAL_FORCE_INTEARACTION_BASE */
