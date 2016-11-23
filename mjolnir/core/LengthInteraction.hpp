#ifndef MJOLNIR_BOND_LENGTH_INTERACTION
#define MJOLNIR_BOND_LENGTH_INTERACTION
#include "Particle.hpp"
#include "PotentialBase.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>

namespace mjolnir
{

template<typename traitsT>
class LengthInteraction
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type> particle_type;
    typedef PotentialBase<traitsT, 2> potential_type;

  public:

    LengthInteraction() = default;
    ~LengthInteraction() = default;

    void
    calc_force(particle_type& p1, particle_type& p2, const potential_type& pot);

    energy_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const potential_type& pot);

};

template<typename traitsT>
inline void
LengthInteraction<traitsT>::calc_force(particle_type& p1, particle_type& p2, 
        const PotentialBase<traitsT, 2>& pot)
{
    const real_type f = -1 * pot.derivative({{p1.position, p2.position}});
    const coordinate_type dpos = p2.position - p1.position;
    const coordinate_type force = dpos * (fast_inv_sqrt(length_sq(dpos)) * f);
    p1.force -= force;
    p2.force += force;
    return;
}

template<typename traitsT>
inline typename LengthInteraction<traitsT>::energy_type
LengthInteraction<traitsT>::calc_energy(
        const particle_type& p1, const particle_type& p2,
        const PotentialBase<traitsT, 2>& pot)
{
    return pot.potential({{p1.position, p2.position}});
}


} // mjolnir
#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
