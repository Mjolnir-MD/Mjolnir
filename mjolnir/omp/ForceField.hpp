#ifndef MJOLNIR_OMP_FORCE_FIELD_HPP
#define MJOLNIR_OMP_FORCE_FIELD_HPP
#include <mjolnir/core/ForceField.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
class ForceField<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type              = OpenMPSimulatorTraits<realT, boundaryT>;
    using real_type                = typename traits_type::real_type;
    using coordinate_type          = typename traits_type::coordinate_type;
    using system_type              = System<traits_type>;
    using local_forcefield_type    = LocalForceField<traits_type>;
    using global_forcefield_type   = GlobalForceField<traits_type>;
    using external_forcefield_type = ExternalForceField<traits_type>;

  public:

    ForceField(local_forcefield_type&& local, global_forcefield_type&& global)
        : local_(std::move(local)), global_(std::move(global))
    {}
    ForceField(local_forcefield_type&&    local,
               global_forcefield_type&&   global,
               external_forcefield_type&& external)
        : local_(std::move(local)), global_(std::move(global)),
          external_(std::move(external))
    {}
    ForceField()  = default;
    ~ForceField() = default;
    ForceField(const ForceField&) = delete;
    ForceField(ForceField&&)      = default;
    ForceField& operator=(const ForceField&) = delete;
    ForceField& operator=(ForceField&&)      = default;

    // this modify system::topology by using local interaction info.
    void initialize(system_type& sys)
    {
        // first, fetch current topology
        local_.write_topology(sys.topology());
        sys.topology().construct_molecules();

        // based on the topology, make exclusion list
        local_   .initialize(sys);
        global_  .initialize(sys);
        external_.initialize(sys);
    }

    // update parameters like temperature, ionic concentration, etc...
    void update(const system_type& sys)
    {
        local_   .update(sys);
        global_  .update(sys);
        external_.update(sys);
    }

    // update margin of neighbor list
    void update_margin(const real_type dmargin, const system_type& sys)
    {
        local_   .update_margin(dmargin, sys);
        global_  .update_margin(dmargin, sys);
        external_.update_margin(dmargin, sys);
    }

    void calc_force(system_type& sys) const noexcept
    {
#pragma omp parallel
        {
            // calc_force uses `nowait` to speedup. To do that, it needs to be
            // inside a parallel region. So here, we need to wrap `calc_force`
            // by a `parallel` region.
            local_   .calc_force(sys);
            global_  .calc_force(sys);
            external_.calc_force(sys);
        }
        // merge thread-local forces into master
        sys.merge_forces();
    }
    real_type calc_energy(const system_type& sys) const noexcept
    {
        return local_.calc_energy(sys) + global_.calc_energy(sys) +
            external_.calc_energy(sys);
    }

    std::string list_energy_name() const
    {
        return local_.list_energy() + global_.list_energy() +
            external_.list_energy();
    }
    std::string dump_energy(const system_type& sys) const
    {
        return local_.dump_energy(sys) + global_.dump_energy(sys) +
            external_.dump_energy(sys);
    }

    local_forcefield_type    const& local()    const noexcept {return local_;}
    local_forcefield_type    &      local()          noexcept {return local_;}
    global_forcefield_type   const& global()   const noexcept {return global_;}
    global_forcefield_type   &      global()         noexcept {return global_;}
    external_forcefield_type const& external() const noexcept {return external_;}
    external_forcefield_type &      external()       noexcept {return external_;}

  private:

    local_forcefield_type    local_;
    global_forcefield_type   global_;
    external_forcefield_type external_;
};

} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
