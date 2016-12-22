#ifndef MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#define MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#include "ParticleContainer.hpp" 
#include "GlobalPotentialBase.hpp"

namespace mjolnir
{

template<typename traitsT>
class GlobalDistanceInteraction
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef GlobalPotentialBase<traits_type> potential_type;

  public:
    GlobalDistanceInteraction() = default;
    ~GlobalDistanceInteraction() = default;

    void
    calc_force(particle_container_type& pcon, potential_type& pot);

    real_type
    calc_energy(const particle_container_type& pcon,
                const potential_type& pot) const;

  private:
    // CellList
};

template<typename traitsT>
inline void
GlobalDistanceInteraction<traitsT>::calc_force(
        particle_container_type& pcon, potential_type& pot)
{
    for(std::size_t i=0; i<pcon.size()-1; ++i)
    {
        for(std::size_t j=i+1; j<pcon.size(); ++j)
        {
            const coordinate_type rij = pcon[j].position - pcon[i].position;
            const real_type       l   = length(rij);
            const coordinate_type f   = rij * (pot.derivative(i, j, l) / l);
            pcon[i].force += f;
            pcon[j].force -= f;
        }
    }
    return ;
}

template<typename traitsT>
inline typename GlobalDistanceInteraction<traitsT>::real_type
GlobalDistanceInteraction<traitsT>::calc_energy(
        const particle_container_type& pcon, const potential_type& pot) const
{
    real_type e = 0.0;
    for(std::size_t i=0; i<pcon.size()-1; ++i)
    {
        for(std::size_t j=i+1; j<pcon.size(); ++j)
        {
            const real_type l = length(pcon[j].position - pcon[i].position);
            e += pot.potential(i, j, l);
        }
    }
    return e;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_DISTANCE_INTEARACTION */
