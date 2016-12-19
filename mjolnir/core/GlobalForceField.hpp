#ifndef MJOLNIR_GLOBAL_FORCE_FIELD
#define MJOLNIR_GLOBAL_FORCE_FIELD
#include "GlobalDistanceInteraction.hpp"
#include "GlobalPotentialBase.hpp"
#include <vector>
#include <array>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class GlobalForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef GlobalPotentialBase<traitsT> potential_type;
    typedef std::unique_ptr<potential_type> potential_ptr;
    typedef GlobalDistanceInteraction<traitsT> distance_interaction_type;

  public:
    GlobalForceField() = default;
    ~GlobalForceField() = default;

    void emplace(potential_ptr&& pot)
    {
        potentials_.emplace_back(std::forward(pot));
    }

    void calc_force(ParticleContainer<traitsT>& pcon);
    real_type calc_force(const ParticleContainer<traitsT>& pcon);

  private:
    distance_interaction_type distance_interaction_;
    std::vector<potential_ptr> potentials_;
};

template<typename traitsT>
void GlobalForceField<traitsT>::calc_force(ParticleContainer<traitsT>& pcon)
{
    for(auto iter = potentials_.cbegin(); iter != potentials_.cend(); ++iter)
        distance_interaction_.calc_force(pcon, *iter);
    return;
}

template<typename traitsT>
typename GlobalForceField<traitsT>::real_type
GlobalForceField<traitsT>::calc_force(const ParticleContainer<traitsT>& pcon)
{
    real_type energy = 0.;
    for(auto iter = potentials_.cbegin(); iter != potentials_.cend(); ++iter)
        energy += distance_interaction_.calc_energy(pcon, *iter);
    return energy;
}

} // mjolnir

#endif /* MJOLNIR_GLOBAL_FORCE_FIELD */
