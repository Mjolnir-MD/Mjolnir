#ifndef MJOLNIR_CORE_EXTERNAL_FORCE_INTEARACTION_BASE_HPP
#define MJOLNIR_CORE_EXTERNAL_FORCE_INTEARACTION_BASE_HPP
#include <mjolnir/core/System.hpp>
#include <string>

namespace mjolnir
{

template<typename traitsT>
class ExternalForceInteractionBase
{
  public:
    using traits_type     = traitsT;
    using system_type     = System<traits_type>;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using boundary_type   = typename traits_type::boundary_type;

  public:
    virtual ~ExternalForceInteractionBase() {}

    virtual void initialize(const system_type&) = 0;
    virtual void update    (const system_type&) = 0;
    virtual void update_margin(const real_type, const system_type&) = 0;

    virtual void      calc_force (system_type&)       const noexcept = 0;
    virtual real_type calc_energy(const system_type&) const noexcept = 0;

    virtual ExternalInteractionBase* clone() const noexcept = 0;

    virtual std::string name() const = 0;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif


} // mjolnir
#endif/* MJOLNIR_EXTERNAL_FORCE_INTEARACTION_BASE */
