#ifndef MJOLNIR_CORE_FORCE_FIELD_BASE_HPP
#define MJOLNIR_CORE_FORCE_FIELD_BASE_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{

// Interface of top-level ForceFields.
template<typename traitsT>
class ForceFieldBase
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;

  public:

    virtual ~ForceFieldBase() = default;

    virtual void initialize(system_type& sys) = 0;
    virtual void update(const system_type& sys) = 0;

    virtual void reduce_margin(const real_type dmargin, const system_type& sys) = 0;
    virtual void  scale_margin(const real_type scale,   const system_type& sys) = 0;

    virtual void calc_force(system_type& sys) const noexcept = 0;
    virtual real_type calc_energy(const system_type& sys) const noexcept = 0;

    virtual std::vector<std::string> list_energy_name() const = 0;
    virtual std::vector<real_type> dump_energy(const system_type& sys) const = 0;
};

} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
