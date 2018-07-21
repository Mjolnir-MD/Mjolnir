#ifndef MJOLNIR_GLOBAL_INTEARACTION_BASE
#define MJOLNIR_GLOBAL_INTEARACTION_BASE
#include <mjolnir/core/System.hpp>

namespace mjolnir
{

template<typename traitsT>
class GlobalInteractionBase
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename traits_type::boundary_type   boundary_type;

  public:

    virtual ~GlobalInteractionBase() = default;

    virtual void initialize(const system_type& sys) = 0;
    virtual void update    (const system_type& sys) = 0;

    virtual void      calc_force (system_type&)             = 0;
    virtual real_type calc_energy(const system_type&) const = 0;

    virtual std::string name() const = 0;
};

} // mjolnir
#endif/* MJOLNIR_GLOBAL_INTEARACTION_BASE */
