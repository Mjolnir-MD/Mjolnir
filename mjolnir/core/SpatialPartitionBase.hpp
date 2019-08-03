#ifndef MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
#define MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>

namespace mjolnir
{
// The following classes would inherits this.
// - NaivePairCalculation
// - VerletList
// - UnlimitedGridCellList
// - PeriodicGridCellList
template<typename traitsT, typename PotentialT>
class SpatialPartitionBase
{
  public:

    using traits_type         = traitsT;
    using system_type         = System<traits_type>;
    using boundary_type       = typename traits_type::boundary_type;
    using real_type           = typename traits_type::real_type;
    using coordinate_type     = typename traits_type::coordinate_type;

    using potential_type      = PotentialT;
    using pair_parameter_type = typename potential_type::pair_parameter_type;
    using neighbor_list_type  = NeighborList<pair_parameter_type>;
    using neighbor_type       = typename neighbor_list_type::neighbor_type;
    using range_type          = typename neighbor_list_type::range_type;

  public:

    SpatialPartitionBase() = default;
    virtual ~SpatialPartitionBase() = default;
    SpatialPartitionBase(SpatialPartitionBase const&) = default;
    SpatialPartitionBase(SpatialPartitionBase&&)      = default;
    SpatialPartitionBase& operator=(SpatialPartitionBase const&) = default;
    SpatialPartitionBase& operator=(SpatialPartitionBase&&)      = default;

    virtual bool valid() const noexcept = 0;

    virtual void initialize(const system_type&, const potential_type&) = 0;

    virtual void make  (const system_type&, const potential_type&) = 0;
    virtual void update(const real_type,    const system_type&,
                        const potential_type&) = 0;

    virtual real_type cutoff() const noexcept = 0;
    virtual real_type margin() const noexcept = 0;

    virtual range_type partners(std::size_t i) const noexcept = 0;
};
} // mjolnir
#endif// MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
