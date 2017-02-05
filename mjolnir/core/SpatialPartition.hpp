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
    typedef std::vector<index_list> verlet_list_type;
    typedef std::vector<index_list> except_list_type;

  public:

    SpatialPartition() noexcept = default;
    virtual ~SpatialPartition() noexcept = default;

    virtual bool valid() const noexcept = 0;

    virtual void make(const particle_container_type& pcon) = 0;
    virtual void update(const particle_container_type& pcon) = 0;
    virtual void update(const particle_container_type& pcon,
                        const time_type dt) = 0;

    void add_except(const index_type i, const index_type j);
    void set_except(const except_list_type& ex);
    void set_except(except_list_type&& ex);

    index_list const& partners(const std::size_t i) const {return list_[i];}
    index_list &      partners(const std::size_t i)       {return list_[i];}

  protected:

    verlet_list_type list_;
    except_list_type except_;
};

template<typename traitsT>
inline void
SpatialPartition<traitsT>::add_except(const index_type i, const index_type j)
{
    this->except_.at(i).push_back(j);
    this->except_.at(j).push_back(i);
    return;
}

template<typename traitsT>
inline void
SpatialPartition<traitsT>::set_except(const except_list_type& ex)
{
    except_ = ex;
    return;
}

template<typename traitsT>
inline void
SpatialPartition<traitsT>::set_except(except_list_type&& ex)
{
    except_ = std::forward<except_list_type>(ex);
    return;
}

} // mjolnir
#endif /*MJOLNIR_CORE_SPATIAL_PARTITIONING*/
