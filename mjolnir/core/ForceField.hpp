#ifndef MJOLNIR_FORCE_FIELD
#define MJOLNIR_FORCE_FIELD
#include "ParticleContainer.hpp"
#include "LocalForceField.hpp"
#include "GlobalForceField.hpp"

namespace mjolnir
{

template<typename traitsT>
class ForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::energy_type energy_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    ForceField(LocalForceField<traitsT>&& local,
               GlobalForceField<traitsT>&& global)
    : local_(std::forward<LocalForceField<traitsT>>(local)),
      global_(std::forward<GlobalForceField<traitsT>>(global))
    {}
    ~ForceField() = default;

    ForceField(const ForceField&) = delete;
    ForceField(ForceField&&)      = default;
    ForceField& operator=(const ForceField&) = delete;
    ForceField& operator=(ForceField&&)      = default;

    void initialize(const ParticleContainer<traitsT>& pcon, const time_type dt);
    void        calc_force(ParticleContainer<traitsT>& pcon);
    energy_type calc_energy(const ParticleContainer<traitsT>& pcon) const;

    void reset_parameter(const std::string& name, const real_type val);

  private:

    LocalForceField<traits_type>  local_;
    GlobalForceField<traits_type> global_;
};
template<typename traitsT>
inline void ForceField<traitsT>::initialize(
        const ParticleContainer<traitsT>& pcon, const time_type dt)
{
//     local_.initialize(pcon, dt);
    global_.initialize(pcon, dt);
    return;
}

template<typename traitsT>
inline void ForceField<traitsT>::calc_force(ParticleContainer<traitsT>& pcon)
{
    this->local_.calc_force(pcon);
    this->global_.calc_force(pcon);
    return ;
}

template<typename traitsT>
inline typename ForceField<traitsT>::energy_type
ForceField<traitsT>::calc_energy(const ParticleContainer<traitsT>& pcon) const
{
    return this->local_.calc_energy(pcon) + this->global_.calc_energy(pcon);
}

template<typename traitsT>
inline void ForceField<traitsT>::reset_parameter(
        const std::string& name, const real_type val)
{
    local_.reset_parameter(name, val);
    global_.reset_parameter(name, val);
    return;
}


} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
