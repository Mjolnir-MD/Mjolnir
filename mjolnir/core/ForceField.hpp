#ifndef MJOLNIR_FORCE_FIELD
#define MJOLNIR_FORCE_FIELD
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/core/GlobalForceField.hpp>
#include <mjolnir/core/ExternalForceField.hpp>

namespace mjolnir
{

template<typename traitsT>
class ForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef System<traits_type>             system_type;
    typedef LocalForceField<traits_type>    local_forcefield_type;
    typedef GlobalForceField<traits_type>   global_forcefield_type;
    typedef ExternalForceField<traits_type> external_forcefield_type;

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
        sys.topology().construct_chains();

        // based on the topology, make exclusion list
           local_.initialize(sys);
          global_.initialize(sys);
        external_.initialize(sys);
    }

    // update parameters like temperature, ionic concentration, etc...
    void update(const system_type& sys)
    {
           local_.update(sys);
          global_.update(sys);
        external_.update(sys);
    }

    void calc_force(system_type& sys)
    {
           local_.calc_force(sys);
          global_.calc_force(sys);
        external_.calc_force(sys);
    }
    real_type calc_energy(const system_type& sys) const
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

  private:

    local_forcefield_type    local_;
    global_forcefield_type   global_;
    external_forcefield_type external_;
};

} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
