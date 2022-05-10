#ifndef MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
#define MJOLNIR_CORE_SPATIAL_PARTITON_BASE_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
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
template<typename traitsT, typename potentialT>
class SpatialPartitionBase
{
  public:

    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using system_type         = System<traits_type>;
    using boundary_type       = typename traits_type::boundary_type;
    using real_type           = typename traits_type::real_type;
    using coordinate_type     = typename traits_type::coordinate_type;

    using parameter_list_type = ParameterListBase<traits_type, potential_type>;
    using neighbor_list_type  = NeighborList<typename potential_type::parameter_type>;
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
            const system_type&, const parameter_list_type&) = 0;

    virtual void make  (neighbor_list_type&,
            const system_type&, const parameter_list_type&) = 0;

    virtual bool reduce_margin(neighbor_list_type&, const real_type,
            const system_type&, const parameter_list_type&) = 0;

    virtual bool scale_margin(neighbor_list_type&, const real_type,
            const system_type&, const parameter_list_type&) = 0;

    virtual real_type cutoff() const noexcept = 0;
    virtual real_type margin() const noexcept = 0;

    virtual SpatialPartitionBase* clone() const = 0;
};

template<typename traitsT, typename potentialT>
class SpatialPartition
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using partition_base_type = SpatialPartitionBase<traits_type, potential_type>;

    using system_type         = typename partition_base_type::system_type        ;
    using boundary_type       = typename partition_base_type::boundary_type      ;
    using real_type           = typename partition_base_type::real_type          ;
    using coordinate_type     = typename partition_base_type::coordinate_type    ;

    using parameter_list_type = typename partition_base_type::parameter_list_type;
    using neighbor_list_type  = typename partition_base_type::neighbor_list_type ;
    using neighbor_type       = typename partition_base_type::neighbor_type      ;
    using range_type          = typename partition_base_type::range_type         ;

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

    void initialize(const system_type& sys, const parameter_list_type& params)
    {
        partition_->initialize(neighbors_, sys, params);
        return;
    }

    void make  (const system_type& sys, const parameter_list_type& params)
    {
        partition_->make(neighbors_, sys, params);
        return ;
    }

    // reduce_margin return true if neighbour list is updated
    bool reduce_margin(const real_type dmargin, const system_type& sys,
                       const parameter_list_type& params)
    {
        return partition_->reduce_margin(neighbors_, dmargin, sys, params);
    }
    bool scale_margin(const real_type scale, const system_type& sys,
                      const parameter_list_type& params)
    {
        return partition_->scale_margin(neighbors_, scale, sys, params);
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
