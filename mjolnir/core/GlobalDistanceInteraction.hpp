#ifndef MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#define MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#include "GlobalInteractionBase.hpp" 
#include "SpatialPartition.hpp" 
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class GlobalDistanceInteraction : public GlobalInteractionBase<traitsT>
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef GlobalPotentialBase<traits_type> potential_type;
    typedef SpatialPartition<traits_type> spartial_partitioning_type;

  public:
    GlobalDistanceInteraction() = default;
    GlobalDistanceInteraction(std::unique_ptr<spatial_partitioning_type>&& sp)
        : spartial_partition_(std::forward<spatial_partitioning_type>(sp))
    {}
    ~GlobalDistanceInteraction() = default;

    void
    calc_force(particle_container_type& pcon, potential_type& pot) override;

    real_type
    calc_energy(const particle_container_type& pcon,
                const potential_type& pot) const override;

    void set_spartial_partition(std::unique_ptr<spatial_partitioning_type>&& sp)
    {
        spartial_partition_ = std::forward<spartial_partitioning_type>(sp);
    }

  private:
    std::unique_ptr<spatial_partitioning_type> spatrial_partition_;
};

template<typename traitsT>
inline void
GlobalDistanceInteraction<traitsT>::calc_force(
        particle_container_type& pcon, potential_type& pot)
{
    spatrial_partition->update(pcon);
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        typename spartial_partitioning_type::index_list const& partners =
            spartial_partition_->partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            // consider periodic boundary case
            const coordinate_type rij = pcon[j].position - pcon[i].position;
            const real_type       l   = length(rij);
            const coordinate_type f   = rij * (pot.derivative(i, j, l) / l);
            pcon[i].force += f;
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
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        typename spartial_partitioning_type::index_list const& partners =
            spartial_partition_->partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = pcon[j].position - pcon[i].position;
            const real_type l = length(rij);
            e += pot.potential(i, j, l);
        }
    }
    return e * 0.5; // double count factor: responsible for spatrial_partition?
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_DISTANCE_INTEARACTION */
