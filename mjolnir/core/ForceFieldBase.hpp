#ifndef MJOLNIR_CORE_FORCE_FIELD_BASE_HPP
#define MJOLNIR_CORE_FORCE_FIELD_BASE_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/core/ConstraintForceField.hpp>

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
    using topology_type   = Topology;
    using constraint_type = ConstraintForceField<traits_type>;

  public:

    virtual ~ForceFieldBase() = default;

    virtual void initialize(const system_type& sys) = 0;
    virtual void update(const system_type& sys) = 0;

    virtual void reduce_margin(const real_type dmargin, const system_type& sys) = 0;
    virtual void  scale_margin(const real_type scale,   const system_type& sys) = 0;

    virtual void calc_force(system_type& sys) const noexcept = 0;
    virtual void calc_force_and_virial(system_type& sys) const noexcept = 0;
    virtual real_type calc_energy(const system_type& sys) const noexcept = 0;

    // format names of all the interactions for .ene file. must not contain '\n'.
    virtual void format_energy_name(std::string&) const = 0;
    // format all energies and related stuff for .ene file. must not contain '\n'.
    // returns the total energy as the return value.
    virtual real_type format_energy(const system_type&, std::string&) const = 0;

    virtual topology_type const& topology() const noexcept = 0;

    // Since constraint should be handled by an integrator in a different way
    // from other potentials, it exposes constraint settings.
    // The reason why constraint is in ForceField is that it affects to the
    // Topology of the forcefield.
    virtual constraint_type const& constraint() const noexcept = 0;
};

} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
