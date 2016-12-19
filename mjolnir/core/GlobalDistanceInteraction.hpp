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
                const potential_type& pot);

  private:
    // CellList
};

template<typename traitsT>
inline void
GlobalDistanceInteraction<traitsT>::calc_force(
        particle_container_type& p1, potential_type& p2)
{
    // TODO
    return ;
}

template<typename traitsT>
inline typename GlobalDistanceInteraction<traitsT>::real_type
GlobalDistanceInteraction<traitsT>::calc_energy(
        const particle_container_type& pcon, const potential_type& pot)
{
    // TODO
    return 0.;
}


} // mjolnir
#endif /* MJOLNIR_GLOBAL_DISTANCE_INTEARACTION */
