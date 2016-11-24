#ifndef MJOLNIR_BOND_LENGTH_INTERACTION
#define MJOLNIR_BOND_LENGTH_INTERACTION
#include "Particle.hpp"
#include "LocalPotentialBase.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>

namespace mjolnir
{

/*! @brief calculate energy and force of Bond length type local interaction */
template<typename traitsT>
class BondLengthInteraction
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type> particle_type;
    typedef LocalPotentialBase<traits_type> potential_type;

  public:

    BondLengthInteraction() = default;
    ~BondLengthInteraction() = default;

    void
    calc_force(particle_type& p1, particle_type& p2, const potential_type& pot);

    real_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const potential_type& pot);
};

template<typename traitsT>
inline void
BondLengthInteraction<traitsT>::calc_force(particle_type& p1, particle_type& p2,
        const potential_type& pot)
{
    const real_type f = -1 * pot.derivative(length(p1.position - p2.position));
    const coordinate_type dpos = p2.position - p1.position;
    const coordinate_type force = dpos * (fast_inv_sqrt(length_sq(dpos)) * f);
    p1.force -= force;
    p2.force += force;
    return;
}

template<typename traitsT>
inline typename BondLengthInteraction<traitsT>::real_type
BondLengthInteraction<traitsT>::calc_energy(
        const particle_type& p1, const particle_type& p2,
        const potential_type& pot)
{
    return pot.potential(length(p1.position - p2.position));
}

} // mjolnir
#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
