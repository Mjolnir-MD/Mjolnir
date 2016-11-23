#ifndef MJOLNIR_GLOBAL_FORCE_FIELD
#define MJOLNIR_GLOBAL_FORCE_FIELD
#include "PotentialBase.hpp"
#include "LengthInteraction.hpp"
#include <utility>
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
    typedef PotentialBase<traitsT, 2> potential_type;
    typedef LengthInteraction<traitsT> length_interaction_type;
    typedef std::unique_ptr<potential_type> potential_ptr;

  public:
    GlobalForceField() = default;
    ~GlobalForceField() = default;

    void emplace(std::size_t i, std::size_t j, potential_ptr&& pot)
    {
        potentials_.emplace_back({{i, j}}, std::forward(pot));
    }

    void calc_force(ParticleContainer& pcon);
    real_type calc_force(const ParticleContainer& pcon);

  private:
    length_interaction_type length_interaction_;
    std::vector<std::pair<std::array<std::size_t, 2>, potential_ptr>> potentials_;
};

template<typename traitsT>
void GlobalForceField<traitsT>::calc_force(ParticleContainer& pcon)
{
    for(auto iter = potentials_.cbegin(); iter != potentials_.cend(); ++iter)
        length_interaction.calc_force(
            pcon.at(iter->first[0]), pcon.at(iter->first[1]), *iter->second);
    return;
}

template<typename traitsT>
real_type GlobalForceField<traitsT>::calc_force(const ParticleContainer& pcon)
{
    real_type energy = 0.;
    for(auto iter = potentials_.cbegin(); iter != potentials_.cend(); ++iter)
        energy += length_interaction.calc_energy(
            pcon.at(iter->first[0]), pcon.at(iter->first[1]), *iter->second);
    return energy;
}

} // mjolnir

#endif /* MJOLNIR_GLOBAL_FORCE_FIELD */
