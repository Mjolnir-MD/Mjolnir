#ifndef MJOLNIR_CORE_SPATIAL_PARTITION_FOR_MU
#define MJOLNIR_CORE_SPATIAL_PARTITION_FOR_MU
#include "System.hpp"
#include <algorithm>
#include <limits>

namespace mjolnir
{
namespace mockup
{
template<typename traitsT>
class SpatialPartitionForMU
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::vector<std::size_t> index_array;
    
  public:

    SpatialPartitionForMU();
    ~SpatialPartitionForMU();
    
};

}//end namspace mockup.
    
}//end namespace mjolnir
#endif// MJOLNIR_CORE_SPATIAL_PARTITION_FOR_MU
