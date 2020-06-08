#ifndef MJOLNIR_CORE_FORCE_FIELD_HPP
#define MJOLNIR_CORE_FORCE_FIELD_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/core/GlobalForceField.hpp>
#include <mjolnir/core/ExternalForceField.hpp>
#include <mjolnir/core/ConstraintForceField.hpp>

namespace mjolnir
{

template<typename traitsT>
class ForceField
{
  public:
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

    ForceField(local_forcefield_type&&    local,
               global_forcefield_type&&   global,
               external_forcefield_type&& external,
               constraint_forcefield_type&& constraint)
        : local_(std::move(local)), global_(std::move(global)),
          external_(std::move(external)), constraint_(std::move(constraint))
    {}

    ForceField()  = default;
    ~ForceField() = default;
    ForceField(const ForceField&) = default;
    ForceField(ForceField&&)      = default;
    ForceField& operator=(const ForceField&) = default;
    ForceField& operator=(ForceField&&)      = default;

    // this modify system::topology by using local interaction info.
    void initialize(system_type& sys)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        MJOLNIR_LOG_INFO("writing current topology");
        topology_.resize(sys.size());
        local_.write_topology(topology_);
        constraint_.write_topology(topology_);
        topology_.construct_molecules();

        MJOLNIR_LOG_INFO("initializing forcefields");
        local_     .initialize(sys);
        global_    .initialize(sys, topology_);
        external_  .initialize(sys);
        return;
    }

    // update parameters like temperature, ionic concentration, etc...
    void update(const system_type& sys)
    {
        local_   .update(sys);
        global_  .update(sys, this->topology_);
        external_.update(sys);
        return;
    }

    // update margin of neighbor list
    void reduce_margin(const real_type dmargin, const system_type& sys)
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
    void scale_margin(const real_type scale, const system_type& sys)
    {
        local_   .scale_margin(scale, sys);
        global_  .scale_margin(scale, sys);
        external_.scale_margin(scale, sys);
        return;
    }

    void calc_force(system_type& sys) const noexcept
    {
        local_   .calc_force(sys);
        global_  .calc_force(sys);
        external_.calc_force(sys);
        return;
    }
    real_type calc_energy(const system_type& sys) const noexcept
    {
        return local_.calc_energy(sys) + global_.calc_energy(sys) +
            external_.calc_energy(sys);
    }

    std::vector<std::string> list_energy_name() const
    {
        auto retval = local_.list_energy();
        auto glo    = global_.list_energy();
        auto ext    = external_.list_energy();

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(std::make_move_iterator(glo.begin()),
                  std::make_move_iterator(glo.end()),
                  std::back_inserter(retval));

        std::copy(std::make_move_iterator(ext.begin()),
                  std::make_move_iterator(ext.end()),
                  std::back_inserter(retval));

        return retval;
    }
    std::vector<real_type> dump_energy(const system_type& sys) const
    {
        auto retval = local_.dump_energy(sys);
        auto glo    = global_.dump_energy(sys);
        auto ext    = external_.dump_energy(sys);

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(glo.begin(), glo.end(), std::back_inserter(retval));
        std::copy(ext.begin(), ext.end(), std::back_inserter(retval));
        return retval;
    }

    topology_type              const& topology()   const noexcept {return topology_;}
    topology_type              &      topology()         noexcept {return topology_;}
    local_forcefield_type      const& local()      const noexcept {return local_;}
    local_forcefield_type      &      local()            noexcept {return local_;}
    global_forcefield_type     const& global()     const noexcept {return global_;}
    global_forcefield_type     &      global()           noexcept {return global_;}
    external_forcefield_type   const& external()   const noexcept {return external_;}
    external_forcefield_type   &      external()         noexcept {return external_;}
    constraint_forcefield_type const& constraint() const noexcept {return constraint_;}
    constraint_forcefield_type &      constraint()       noexcept {return constraint_;}

  private:

    topology_type               topology_;
    local_forcefield_type       local_;
    global_forcefield_type      global_;
    external_forcefield_type    external_;
    constraint_forcefield_type  constraint_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own specialization to avoid data race.
    // So this implementation should not be instanciated with the OpenMP traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif

};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ForceField<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class ForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
