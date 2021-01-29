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

    virtual bool reduce_margin(neighbor_list_type&, const real_type,
            const system_type&, const potential_type&) = 0;

    virtual void scale_margin(neighbor_list_type&, const real_type,
            const system_type&, const potential_type&) = 0;

    virtual real_type cutoff() const noexcept = 0;
    virtual real_type margin() const noexcept = 0;

    virtual SpatialPartitionBase* clone() const = 0;
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

    using partition_base_type = SpatialPartitionBase<traits_type, potential_type>;
    using partition_type      = std::unique_ptr<partition_base_type>;

  public:

    explicit SpatialPartition(partition_type&& part)
        : partition_(std::move(part))
    {}
    ~SpatialPartition() = default;
    SpatialPartition(SpatialPartition&&)            = default;
    SpatialPartition& operator=(SpatialPartition&&) = default;

    SpatialPartition(SpatialPartition const& other)
        : partition_(other.base().clone()), neighbors_(other.neighbors())
    {}
    SpatialPartition& operator=(SpatialPartition const& other)
    {
        this->partition_.reset(other.base().clone());
        this->neighbors_ = other.neighbors();
        return *this;
    }

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

    // reduce_margin return true if neighbour list is updated
    bool reduce_margin(const real_type dmargin, const system_type& sys,
                       const potential_type& pot)
    {
        return partition_->reduce_margin(neighbors_, dmargin, sys, pot);
    }
    void scale_margin(const real_type scale, const system_type& sys,
                      const potential_type& pot)
    {
        partition_->scale_margin(neighbors_, scale, sys, pot);
        return ;
    }

    real_type cutoff() const noexcept {return partition_->cutoff();}
    real_type margin() const noexcept {return partition_->margin();}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

    // for testing
    partition_base_type const& base() const noexcept {return *partition_;}
    partition_base_type &      base()       noexcept {return *partition_;}

    neighbor_list_type const& neighbors() const noexcept {return neighbors_;}
    neighbor_list_type &      neighbors()       noexcept {return neighbors_;}

  private:

    partition_type     partition_;
    neighbor_list_type neighbors_;
};

} // mjolnir
#endif// MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
