#ifndef MJOLNIR_NEIGHBOR_LIST_HPP
#define MJOLNIR_NEIGHBOR_LIST_HPP
#include <mjolnir/util/range.hpp>
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

namespace mjolnir
{

//! a container for verlet-list, cell-list, and other spatial indexing methods.
//! the objectives are
//! - store indices of possible partners that has interaction.
//! - store pre-calculated parameters. e.g. epsilon_ij = sqrt(eps_i * eps_j)
//!   for Lennard-Jones takes relatively large computational cost to obtain.
//!   calculate it for all the possible pairs to make force calculation fast.
template<typename parameterT>
class NeighborList
{
  public:
    using parameter_type = parameterT;
    using neighbor_type  = std::pair<std::size_t, parameter_type>;
    using container_type = std::vector<neighbor_type>;
    using range_type     = range<typename container_type::const_iterator>;

  public:

    NeighborList() = default;
    ~NeighborList() = default;
    NeighborList(const NeighborList&) = default;
    NeighborList(NeighborList&&)      = default;
    NeighborList& operator=(const NeighborList&) = default;
    NeighborList& operator=(NeighborList&&)      = default;

    void clear()
    {
        this->neighbors_.clear();
        this->ranges_.clear();
        return;
    }

    void reserve(std::size_t Nparticle, std::size_t Nneighbor)
    {
        this->neighbors_.reserve(Nparticle * Nneighbor);
        this->ranges_   .reserve(Nparticle);
        return;
    }

    template<typename Iterator>
    void add_list_for(const std::size_t i, Iterator first, Iterator last)
    {
        if(this->ranges_.size() <= i)
        {
            this->ranges_.resize(i+1, {0,0});
        }
        ranges_[i].first  = neighbors_.size();
        ranges_[i].second = neighbors_.size() + std::distance(first, last);

        // this might break iterator.
        // this->ranges_ cannot contain pair of iteraor but indices.
        std::copy(first, last, std::back_inserter(this->neighbors_));
        return ;
    }

    range_type operator[](const std::size_t i) const noexcept
    {
        return range_type{
            this->neighbors_.begin() + this->ranges_[i].first,
            this->neighbors_.begin() + this->ranges_[i].second
        };
    }
    range_type at(const std::size_t i) const
    {
        return range_type{
            this->neighbors_.begin() + this->ranges_.at(i).first,
            this->neighbors_.begin() + this->ranges_.at(i).second
        };
    }

  private:
    container_type neighbors_;
    std::vector<std::pair<std::size_t, std::size_t>> ranges_;
    // consider to use flat_map when a few particle interacts?
};

// TODO: when no parameter is needed...

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
