#ifndef MJOLNIR_CORE_GLOBAL_INTEARACTION_BASE_HPP
#define MJOLNIR_CORE_GLOBAL_INTEARACTION_BASE_HPP
#include <mjolnir/core/System.hpp>
#include <string>

namespace mjolnir
{

template<typename traitsT>
class GlobalInteractionBase
{
  public:
    using traits_type     = traitsT;
    using system_type     = System<traits_type>;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using boundary_type   = typename traits_type::boundary_type;
    using topology_type   = Topology;

  public:
    virtual ~GlobalInteractionBase() {}

    virtual void initialize(const system_type&, const topology_type&) = 0;
    virtual void update    (const system_type&, const topology_type&) = 0;
    virtual void reduce_margin(const real_type, const system_type&) = 0;
    virtual void  scale_margin(const real_type, const system_type&) = 0;

    virtual void      calc_force           (system_type&)    const noexcept = 0;
    virtual void      calc_force_and_virial(system_type&)    const noexcept = 0;
    virtual real_type calc_energy(const system_type&)        const noexcept = 0;
    virtual real_type calc_force_and_energy(system_type&)    const noexcept = 0;
    virtual real_type calc_force_virial_energy(system_type&) const noexcept = 0;

    virtual GlobalInteractionBase* clone() const = 0;

    virtual std::string name() const = 0;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif/* MJOLNIR_GLOBAL_INTEARACTION_BASE */
