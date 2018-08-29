#ifndef MJOLNIR_NEIGHBOR_LIST_HPP
#define MJOLNIR_NEIGHBOR_LIST_HPP
#include <mjolnir/util/range.hpp>
#include <mjolnir/util/empty.hpp>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

//! a container for verlet-list, cell-list, and other spatial indexing methods.
//! the objectives are
//! - store indices of possible partners that has interaction.
//! - store pre-calculated parameters. e.g. epsilon_ij = sqrt(eps_i * eps_j)
//!   for Lennard-Jones takes relatively large computational cost to obtain.
//!   calculate it for all the possible pairs to make force calculation fast.

namespace mjolnir
{

// google with "empty base optimization(EBO)".
// by using EBO, it compress the object size when `paramT` is a empty class.
// without this, empty parameter type consumes 8 byte (on x86_64) because
// std::size_t requires 8-byte alignment.

namespace detail
{
// In order to handle `double` or other non-class object, neighbor_element need
// another layer that checks `std::is_empty`. If empty, it inherits it and EBO
// removes the overhead. If not empty(like `double`), it keeps storage for the
// value. `neighbor_elmenet_impl` stores value if `std::true_type` is passed.
template<typename paramT, bool IsEmpty>
struct neighbor_element_impl;

template<typename paramT>
struct neighbor_element_impl<paramT, true> : private paramT // for EBO
{
    using parameter_type = paramT;

    neighbor_element_impl(const parameter_type& p) {}
    neighbor_element_impl(parameter_type&& p)      {}

    parameter_type&       parameter()       noexcept {return *this;}
    parameter_type const& parameter() const noexcept {return *this;}
};

template<typename paramT>
struct neighbor_element_impl<paramT, false>
{
    using parameter_type = paramT;

    neighbor_element_impl(const parameter_type& p) : param(p) {}
    neighbor_element_impl(parameter_type&& p)      : param(std::move(p)) {}

    parameter_type&       parameter()       noexcept {return param;}
    parameter_type const& parameter() const noexcept {return param;}

    paramT param;
};

} // detail

template<typename paramT>
struct neighbor_element
: private detail::neighbor_element_impl<paramT, std::is_empty<paramT>::value>
{
    using base_type =
        detail::neighbor_element_impl<paramT, std::is_empty<paramT>::value>;
    using parameter_type = typename base_type::parameter_type;

    neighbor_element(std::size_t idx, const paramT& p)
        : base_type(p), index(idx)
    {}
    neighbor_element(std::size_t idx, paramT&& p)
        : base_type(std::move(p)), index(idx)
    {}

    neighbor_element()  = default;
    ~neighbor_element() = default;
    neighbor_element(const neighbor_element&) = default;
    neighbor_element(neighbor_element&&)      = default;
    neighbor_element& operator=(const neighbor_element&) = default;
    neighbor_element& operator=(neighbor_element&&)      = default;

    parameter_type&       parameter()       noexcept {return base_type::parameter();}
    parameter_type const& parameter() const noexcept {return base_type::parameter();}

    std::size_t index;
};

static_assert(sizeof(std::size_t) == sizeof(neighbor_element<empty_t>));

template<typename paramT>
inline bool operator==(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return lhs.index == rhs.index;
}
template<typename paramT>
inline bool operator!=(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return lhs.index != rhs.index;
}
template<typename paramT>
inline bool operator<(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return lhs.index < rhs.index;
}
template<typename paramT>
inline bool operator>(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return lhs.index > rhs.index;
}
template<typename paramT>
inline bool operator<=(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return lhs.index <= rhs.index;
}
template<typename paramT>
inline bool operator>=(
    const neighbor_element<paramT>& lhs, const neighbor_element<paramT>& rhs)
{
    return lhs.index >= rhs.index;
}

template<typename parameterT>
class NeighborList
{
  public:
    using parameter_type = parameterT;
    using neighbor_type  = neighbor_element<parameter_type>;
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

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
