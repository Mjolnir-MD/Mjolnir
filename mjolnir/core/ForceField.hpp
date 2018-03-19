#ifndef MJOLNIR_FORCE_FIELD
#define MJOLNIR_FORCE_FIELD
#include "LocalForceField.hpp"
#include "GlobalForceField.hpp"
#include "ExternalForceField.hpp"

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

    ForceField(local_forcefield_type&& local, global_forcefield_type&& global,
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

    void initialize(const system_type& sys, const real_type dt)
    {
        global_.initialize(sys, dt);
    }

    // update parameters like delta t, temperature, ionic concentration, etc...
    void update(const system_type& sys, const real_type dt)
    {
           local_.update(sys, dt);
          global_.update(sys, dt);
        external_.update(sys, dt);
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
