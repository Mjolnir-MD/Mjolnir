#ifndef MJOLNIR_CORE_MULTIPLE_BASIN_UNIT_BASE_HPP
#define MJOLNIR_CORE_MULTIPLE_BASIN_UNIT_BASE_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/string.hpp>
#include <algorithm>
#include <numeric>
#include <memory>

namespace mjolnir
{

// MultipleBasinForceField
//
// In some cases, protein has several domains that undergoes conformational
// changes independently.
//
//  .-.  .-.      __  .-.
// ( A )( B ) -> (A'|( B )
//  `-'  `-'      `-  `-'
//      |            |
//      v            v
//  .-.  __        __ __
// ( A )|B')  ->  (A'|B')
//  `-'  -'        `- -'
//
// We call those A and B as a "unit" of MultipleBasinForceField.
// This class defines the interface of a unit.
//
template<typename traitsT>
class MultipleBasinUnitBase
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using topology_type   = Topology;

  public:

    virtual ~MultipleBasinUnitBase() = default;

    virtual void write_topology(const system_type&, topology_type&) const = 0;
    virtual void initialize    (const system_type&, const topology_type&) = 0;

    virtual void      calc_force (      system_type&) const noexcept = 0;
    virtual real_type calc_energy(const system_type&) const noexcept = 0;

    // update system parameters (e.g. temperature, ionic_strength, etc.)
    virtual void update(const system_type&, const topology_type&) = 0;

    // update margin of neighbor list
    virtual void reduce_margin(const real_type dmargin, const system_type&) = 0;
    virtual void scale_margin (const real_type scale,   const system_type&) = 0;

    // energy output
    virtual std::vector<std::string> list_energy_name() const = 0;
    virtual std::vector<real_type>   dump_energy(const system_type&) const = 0;
};

} // mjolnir
#endif// MJOLNIR_CORE_MULTIPLE_BASIN_FORCE_FIELD_HPP
