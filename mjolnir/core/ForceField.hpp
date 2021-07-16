#ifndef MJOLNIR_CORE_FORCE_FIELD_HPP
#define MJOLNIR_CORE_FORCE_FIELD_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/core/GlobalForceField.hpp>
#include <mjolnir/core/ExternalForceField.hpp>
#include <mjolnir/core/ConstraintForceField.hpp>

namespace mjolnir
{

template<typename traitsT>
class ForceField final : public ForceFieldBase<traitsT>
{
  public:

    using base_type                   = ForceFieldBase<traitsT>;
    using traits_type                 = traitsT;
    using real_type                   = typename traits_type::real_type;
    using coordinate_type             = typename traits_type::coordinate_type;
    using system_type                 = System<traits_type>;
    using topology_type               = Topology;
    using local_forcefield_type       = LocalForceField<traits_type>;
    using global_forcefield_type      = GlobalForceField<traits_type>;
    using external_forcefield_type    = ExternalForceField<traits_type>;
    using constraint_forcefield_type  = ConstraintForceField<traits_type>;

  public:

    ForceField(local_forcefield_type&&      local,
               global_forcefield_type&&     global,
               external_forcefield_type&&   external,
               constraint_forcefield_type&& constraint)
        : local_(std::move(local)), global_(std::move(global)),
          external_(std::move(external)), constraint_(std::move(constraint))
    {}

    ForceField()  = default;
    ~ForceField() override = default;
    ForceField(const ForceField&) = default;
    ForceField(ForceField&&)      = default;
    ForceField& operator=(const ForceField&) = default;
    ForceField& operator=(ForceField&&)      = default;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        MJOLNIR_LOG_INFO("writing current topology");
        topology_.resize(sys.size());
        local_     .write_topology(topology_);
        constraint_.write_topology(topology_);
        topology_.construct_molecules();

        MJOLNIR_LOG_INFO("initializing forcefields");
        local_     .initialize(sys);
        global_    .initialize(sys, topology_);
        external_  .initialize(sys);
        return;
    }

    // update parameters like temperature, ionic concentration, etc...
    void update(const system_type& sys) override
    {
        local_   .update(sys);
        global_  .update(sys, this->topology_);
        external_.update(sys);
        return;
    }

    // update margin of neighbor list
    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        // TODO:
        // dmargin is a maximum (safe) length to reduce the margin. Since global
        // forcefields calculate forces between particles, the safe length is
        // mostly 2 * (largest displacement in the step). It is dmargin.
        // However, since external forcefields calculate forces between particle
        // and an external field, the safe length can be just the largest
        // displacement in the step. The current implementation overestimate the
        // safe length for the external forcefields (it's correct for global
        // forcefields). In most cases, the most time-consuming part is global
        // forcefields, so I implemented in this way for now. If some idea that
        // works more efficiently is came up, this part would be re-implemented.
        local_   .reduce_margin(dmargin, sys);
        global_  .reduce_margin(dmargin, sys);
        external_.reduce_margin(dmargin, sys);
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        local_   .scale_margin(scale, sys);
        global_  .scale_margin(scale, sys);
        external_.scale_margin(scale, sys);
        return;
    }

    void calc_force(system_type& sys) const noexcept override
    {
        sys.preprocess_forces();
        local_   .calc_force(sys);
        global_  .calc_force(sys);
        external_.calc_force(sys);
        sys.postprocess_forces();
        return;
    }
    real_type calc_energy(const system_type& sys) const noexcept override
    {
        return local_.calc_energy(sys) + global_.calc_energy(sys) +
            external_.calc_energy(sys);
    }

    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        real_type energy(0);
        sys.preprocess_forces();
        energy += local_   .calc_force_and_energy(sys);
        energy += global_  .calc_force_and_energy(sys);
        energy += external_.calc_force_and_energy(sys);
        sys.postprocess_forces();
        return energy;
    }

    void format_energy_name(std::string& fmt) const override
    {
        local_   .format_energy_name(fmt);
        global_  .format_energy_name(fmt);
        external_.format_energy_name(fmt);
        return ;
    }
    real_type format_energy(const system_type& sys, std::string& fmt) const override
    {
        real_type total = 0.0;
        total += local_   .format_energy(sys, fmt);
        total += global_  .format_energy(sys, fmt);
        total += external_.format_energy(sys, fmt);
        return total;
    }

    topology_type const& topology() const noexcept override {return topology_;}

    local_forcefield_type      const& local()      const noexcept {return local_;}
    local_forcefield_type      &      local()            noexcept {return local_;}
    global_forcefield_type     const& global()     const noexcept {return global_;}
    global_forcefield_type     &      global()           noexcept {return global_;}
    external_forcefield_type   const& external()   const noexcept {return external_;}
    external_forcefield_type   &      external()         noexcept {return external_;}
    constraint_forcefield_type const& constraint() const noexcept override {return constraint_;}
    constraint_forcefield_type &      constraint()       noexcept {return constraint_;}

  private:

    topology_type               topology_;
    local_forcefield_type       local_;
    global_forcefield_type      global_;
    external_forcefield_type    external_;
    constraint_forcefield_type  constraint_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ForceField<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class ForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
