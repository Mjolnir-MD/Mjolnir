#ifndef MJOLNIR_PAIR_DISTANCE_INTERACTION
#define MJOLNIR_PAIR_DISTANCE_INTERACTION
#include "Particle.hpp"
#include "GlobalPotentialBase.hpp"
#include <cmath>

namespace mjolnir
{

template<typename traitsT, typename valueT>
class PairDistanceInteraction
{
  public:
    typedef traitsT traits_type;
    typedef valueT value_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type> particle_type;
    typedef GlobalPotentialBase<traits_type, value_type> potential_type;

  public:
    PairDistanceInteraction() = default;
    ~PairDistanceInteraction() = default;

    void
    calc_force(particle_type& p1, particle_type& p2,
               const value_type& v1, const value_type& v2,
               const potential_type& pot);

    real_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const value_type& v1, const value_type& v2,
                const potential_type& pot);
};

template<typename traitsT, typename valueT>
void PairDistanceInteraction<traitsT, valueT>::calc_force(
        particle_type& p1, particle_type& p2,
        const value_type& v1, const value_type& v2, const potential_type& pot)
{
    const coordinate_type dpos = p2.position - p1.position;
    const real_type lensq = length_sq(dpos);
    const real_type f = -1 * pot.derivative(v1, v2, std::sqrt(lensq));
    const coordinate_type force = dpos * (fast_inv_sqrt(lensq) * f);
    p1.force -= force;
    p2.force += force;
    return;
}

template<typename traitsT, typename valueT>
typename PairDistanceInteraction<traitsT, valueT>::real_type
PairDistanceInteraction<traitsT, valueT>::calc_energy(
        const particle_type& p1, const particle_type& p2,
        const value_type& v1, const value_type& v2, const potential_type& pot)
{
    return pot.potential(v1, v2, length(p1.position - p2.position));
}

} // mjolnir
#endif /* MJOLNIR_PAIR_DISTANCE_INTERACTION */
