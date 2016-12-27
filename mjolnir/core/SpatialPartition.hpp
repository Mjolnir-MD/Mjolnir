#ifndef MJOLNIR_CORE_SPATIAL_PARTITIONING
#define MJOLNIR_CORE_SPATIAL_PARTITIONING
#include "ParticleContainer.hpp"

namespace mjolnir
{

template<typename traitsT>
class SpatialPartition
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef typename particle_container_type::index_type index_type;
    typedef std::vector<index_type> index_list;
 
  public:

    SpatialPartition()  noexcept = default;
    virtual ~SpatialPartition() noexcept = default;

    virtual bool valid() const noexcept = 0;
    virtual void make(const particle_container_type pcon) = 0;
    virtual void update(const particle_container_type& pcon) = 0;
    virtual void update(const particle_container_type& pcon, const time_type dt) = 0;

    virtual index_list &      partners(const index_type i)       = 0;
    virtual index_list const& partners(const index_type i) const = 0;
};

} // mjolnir
#endif /*MJOLNIR_CORE_SPATIAL_PARTITIONING*/
