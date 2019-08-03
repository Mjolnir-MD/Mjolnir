#ifndef MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
#define MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <memory>

namespace mjolnir
{

// SpatialPartitionBase is a base class of the following classes.
// - NaivePairCalculation
// - VerletList
// - UnlimitedGridCellList
// - PeriodicGridCellList
// It constructs NeighborLists. GlobalInteraction uses this internally.
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

    virtual void initialize(neighbor_list_type&,
            const system_type&, const potential_type&) = 0;

    virtual void make  (neighbor_list_type&,
            const system_type&, const potential_type&) = 0;

    virtual void update(neighbor_list_type&,
            const real_type,    const system_type&, const potential_type&) = 0;

    virtual real_type cutoff() const noexcept = 0;
    virtual real_type margin() const noexcept = 0;
};

template<typename traitsT, typename PotentialT>
class SpatialPartition
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

    using partition_type      =
        std::unique_ptr<SpatialPartitionBase<traits_type, potential_type>>;

  public:

    SpatialPartition(partition_type&& part): partition_(std::move(part)){}

    ~SpatialPartition() = default;
    SpatialPartition(SpatialPartition const&) = default;
    SpatialPartition(SpatialPartition&&)      = default;
    SpatialPartition& operator=(SpatialPartition const&) = default;
    SpatialPartition& operator=(SpatialPartition&&)      = default;

    bool valid() const noexcept {return partition_->valid();}

    void initialize(const system_type& sys, const potential_type& pot)
    {
        partition_->initialize(neighbors_, sys, pot);
        return;
    }

    void make  (const system_type& sys, const potential_type& pot)
    {
        partition_->make(neighbors_, sys, pot);
        return ;
    }
    void update(const real_type dmargin, const system_type& sys,
                const potential_type& pot)
    {
        partition_->update(neighbors_, dmargin, sys, pot);
        return ;
    }

    real_type cutoff() const noexcept {return partition_->cutoff();}
    real_type margin() const noexcept {return partition_->margin();}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    partition_type     partition_;
    neighbor_list_type neighbors_;
};



} // mjolnir
#endif// MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
