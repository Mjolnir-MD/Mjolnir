#ifndef MJOLNIR_LOCAL_INTEARACTION_BASE
#define MJOLNIR_LOCAL_INTEARACTION_BASE
#include <mjolnir/core/System.hpp>
#include <string>

namespace mjolnir
{

template<typename traitsT>
class LocalInteractionBase
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename system_type::topology_type   topology_type;
    typedef typename topology_type::connection_kind_type connection_kind_type;

  public:

    virtual ~LocalInteractionBase() = default;

    virtual void initialize(const system_type&) = 0;
    virtual void update    (const system_type&) = 0;
    virtual void update_margin(const real_type, const system_type&) = 0;

    virtual void      calc_force (system_type&)       const noexcept = 0;
    virtual real_type calc_energy(const system_type&) const noexcept = 0;

    virtual std::string name() const = 0;
    virtual void write_topology(topology_type&) const = 0;
};

} // mjolnir
#endif/* MJOLNIR_GLOBAL_INTEARACTION_BASE */
