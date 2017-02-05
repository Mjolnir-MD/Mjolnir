#ifndef MJOLNIR_INTER_CHAIN_DISTANCE_INTEARACTION
#define MJOLNIR_INTER_CHAIN_DISTANCE_INTEARACTION
#include "GlobalInteractionBase.hpp"
#include "SpatialPartition.hpp"
#include "BoundaryCondition.hpp"
#include <memory>

namespace mjolnir
{

/*! @brief only for inter-chain interaction */
template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class InterChainDistanceInteraction : public GlobalInteractionBase<traitsT>
{
  public:

    typedef traitsT traits_type;
    typedef boundaryT boundary_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef GlobalPotentialBase<traits_type> potential_type;
    typedef SpatialPartition<traits_type> spatial_partitioning_type;
    typedef std::size_t index_type;
    typedef std::vector<index_type> index_list_type;

  public:
    InterChainDistanceInteraction() = default;
    ~InterChainDistanceInteraction() = default;

    InterChainDistanceInteraction(std::unique_ptr<spatial_partitioning_type>&& sp)
        : spatial_partition_(
                std::forward<std::unique_ptr<spatial_partitioning_type>>(sp))
    {}

    InterChainDistanceInteraction(std::unique_ptr<spatial_partitioning_type>&& sp,
                              const coordinate_type& system_size)
        : spatial_partition_(
                std::forward<std::unique_ptr<spatial_partitioning_type>>(sp)),
          boundary(system_size)
    {}

    void add_chain(index_list_type&&      chain);
    void add_chain(const index_list_type& chain);

    void
    calc_force(particle_container_type& pcon, potential_type& pot) override;

    real_type
    calc_energy(const particle_container_type& pcon,
                const potential_type& pot) const override;

    void set_spatial_partition(std::unique_ptr<spatial_partitioning_type>&& sp)
    {
        spatial_partition_ = std::forward<spatial_partitioning_type>(sp);
    }

    void initialize(const particle_container_type& pcon, const time_type dt) override;

  private:

    std::vector<index_list_type> chains_;

    std::unique_ptr<spatial_partitioning_type> spatial_partition_;
    boundary_type boundary;
};

template<typename traitsT, typename boundaryT>
inline void
InterChainDistanceInteraction<traitsT, boundaryT>::add_chain(
        index_list_type&& chain)
{
    chains_.emplace_back(std::forward<index_list_type>(chain));
    return;
}

template<typename traitsT, typename boundaryT>
inline void
InterChainDistanceInteraction<traitsT, boundaryT>::add_chain(
        const index_list_type& chain)
{
    chains_.push_back(chain);
    return ;
}

template<typename traitsT, typename boundaryT>
void InterChainDistanceInteraction<traitsT, boundaryT>::calc_force(
        particle_container_type& pcon, potential_type& pot)
{
    spatial_partition_->update(pcon);
    for(auto iter = chains_.cbegin(); iter != chains_.cend(); ++iter)
    {
        for(auto index = iter->cbegin(); index != iter->cend(); ++index)
        {
            const std::size_t i = *index;
            typename spatial_partitioning_type::index_list const& partners =
                spatial_partition_->partners(i);
            for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
            {
                const std::size_t j = *iter;
                // XXX: for short-range force
                const coordinate_type rij =
                    boundary(pcon[j].position - pcon[i].position);
                const real_type       l = length(rij);
                const coordinate_type f = rij * (pot.derivative(i, j, l) / l);
                pcon[i].force += f;
                pcon[j].force -= f;
            }
        }
    }
    return ;
}

template<typename traitsT, typename boundaryT>
typename InterChainDistanceInteraction<traitsT, boundaryT>::real_type
InterChainDistanceInteraction<traitsT, boundaryT>::calc_energy(
        const particle_container_type& pcon, const potential_type& pot) const
{
    real_type e = 0.0;
    for(auto iter = chains_.cbegin(); iter != chains_.cend(); ++iter)
    {
        for(auto index = iter->cbegin(); index != iter->cend(); ++index)
        {
            const std::size_t i = *index;
            typename spatial_partitioning_type::index_list const& partners =
                spatial_partition_->partners(i);
            for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
            {
                const std::size_t j = *iter;
                // XXX: for short-range force
                const coordinate_type rij =
                    boundary(pcon[j].position - pcon[i].position);
                const real_type l = length(rij);
                e += pot.potential(i, j, l);
            }
        }
    }
    return e;
}

template<typename traitsT, typename boundaryT>
inline void InterChainDistanceInteraction<traitsT, boundaryT>::initialize(
        const particle_container_type& pcon, const time_type dt)
{
    this->spatial_partition_->update(pcon, dt);
    return ;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_DISTANCE_INTEARACTION */
