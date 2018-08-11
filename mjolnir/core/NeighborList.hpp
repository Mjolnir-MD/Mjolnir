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
class NeighborList
{
  public:
    typedef range<typename std::vector<std::size_t>::const_iterator> range_type;

  public:

    NeighborList() = default;
    ~NeighborList() = default;
    NeighborList(const NeighborList&) = default;
    NeighborList(NeighborList&&)      = default;
    NeighborList& operator=(const NeighborList&) = default;
    NeighborList& operator=(NeighborList&&)      = default;

    void clear() {this->idxs_.clear(); this->ranges_.clear(); return;}

    void reserve(std::size_t Nparticle, std::size_t Nneighbor)
    {
        this->idxs_  .reserve(Nparticle * Nneighbor);
        this->ranges_.reserve(Nparticle);
    }

    template<typename Iterator>
    void add_list_for(const std::size_t i, Iterator first, Iterator last)
    {
        if(this->ranges_.size() <= i)
        {
            this->ranges_.resize(i+1, {0,0});
        }
        ranges_[i].first  = idxs_.size();
        ranges_[i].second = idxs_.size() + std::distance(first, last);

        // this might break iterator.
        // this->ranges_ cannot contain pair of iteraor but indices.
        std::copy(first, last, std::back_inserter(this->idxs_));
        return ;
    }

    range_type operator[](const std::size_t i) const noexcept
    {
        return range_type{
            this->idxs_.begin() + this->ranges_[i].first,
            this->idxs_.begin() + this->ranges_[i].second
        };
    }
    range_type at(const std::size_t i) const
    {
        return range_type{
            this->idxs_.begin() + this->ranges_.at(i).first,
            this->idxs_.begin() + this->ranges_.at(i).second
        };
    }

  private:
    std::vector<std::size_t> idxs_;
    std::vector<std::pair<std::size_t, std::size_t>> ranges_;
};



} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
