#ifndef MJOLNIR_CORE_NEIGHBOR_LIST_HPP
#define MJOLNIR_CORE_NEIGHBOR_LIST_HPP
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
// By using EBO, we can compress the object size when `paramT` is a empty class.
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

    explicit neighbor_element_impl(const parameter_type&) {}
    explicit neighbor_element_impl(parameter_type&&)      {}

    neighbor_element_impl()  = default;
    ~neighbor_element_impl() = default;
    neighbor_element_impl(const neighbor_element_impl&) = default;
    neighbor_element_impl(neighbor_element_impl&&)      = default;
    neighbor_element_impl& operator=(const neighbor_element_impl&) = default;
    neighbor_element_impl& operator=(neighbor_element_impl&&)      = default;

    parameter_type&       parameter()       noexcept {return *this;}
    parameter_type const& parameter() const noexcept {return *this;}
};

template<typename paramT>
struct neighbor_element_impl<paramT, false>
{
    using parameter_type = paramT;

    explicit neighbor_element_impl(const parameter_type& p): param(p) {}
    explicit neighbor_element_impl(parameter_type&& p): param(std::move(p)) {}

    neighbor_element_impl()  = default;
    ~neighbor_element_impl() = default;
    neighbor_element_impl(const neighbor_element_impl&) = default;
    neighbor_element_impl(neighbor_element_impl&&)      = default;
    neighbor_element_impl& operator=(const neighbor_element_impl&) = default;
    neighbor_element_impl& operator=(neighbor_element_impl&&)      = default;

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

    // paramT param; // derived from neighbor_element_impl
    std::size_t index;
};

// Check the EBO works and the size of neighbor_element with empty class
// is equal to the size of `std::size_t index;`.
static_assert(sizeof(std::size_t) == sizeof(neighbor_element<empty_t>),
              "checking neighbor_element reduces size of empty object");

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

// ----------------------------------------------------------------------------
// neighbor list

template<typename parameterT>
class NeighborList
{
  public:
    using parameter_type = parameterT;
    using neighbor_type  = neighbor_element<parameter_type>;
    using container_type = std::vector<neighbor_type>;
    using range_type     = range<typename container_type::const_iterator>;

    // XXX this contains the list in the following way.
    //
    // 1. The neighboring list that has interaction partners of each particle
    //    and the parameter for each pairs (e.g. `q_i*q_j` for electrostatics).
    //
    //  __partner of p1__ __partner of p2___________ ____ ... pN__
    // '                 '                          '             '
    // |{p2,q12}|{p3,q13}|{p4,q24}|{p6,q26}|{p7,q27}|{    ...    }|
    //   0th      1st      2nd      3rd      4th      5th ... M-1 th element
    //
    // 2. The range list that contains from where to where the partners are
    //    contained.
    //  +---- partner of p1 starts from 0th element of the above array
    //  | +-- partner of p2 starts from 2nd element of the above array
    //  v v
    // |0|2|5|...|M|
    //    ^ ^     ^
    //    +-+-----+-- the last partner of p1 is 2-1 = 1.
    //      +-----+-- the last pertner of p2 is 5-1 = 4.
    //            +-- the last element of pN is M-1.
    //
    // This list should have N+1 elements because i-th and (i+1)-th elements of
    // the list is needed to obtain the partner.
    //
    // By doing this, we can reduce the memory resource to have the list.

  public:

    NeighborList()  = default;
    ~NeighborList() = default;
    NeighborList(const NeighborList&) = default;
    NeighborList(NeighborList&&)      = default;
    NeighborList& operator=(const NeighborList&) = default;
    NeighborList& operator=(NeighborList&&)      = default;

    void clear()
    {
        this->neighbors_.clear();
        this->ranges_   .clear();
        return;
    }

    void reserve(std::size_t Nparticle, std::size_t Nneighbor)
    {
        this->neighbors_.reserve(Nparticle * Nneighbor);
        this->ranges_   .reserve(Nparticle + 1);
        return;
    }

    // assign the partners in range [first, last) as the partner of particle i.
    template<typename Iterator>
    void add_list_for(const std::size_t i, Iterator first, Iterator last)
    {
        static_assert(std::is_same<
            typename std::iterator_traits<Iterator>::value_type, neighbor_type
            >::value, "The iterator must points neighbor type.");

        if(this->ranges_.size() <= i+1)
        {
            this->ranges_.resize(i+2, 0);
        }
        this->ranges_[i]   = this->neighbors_.size();
        this->ranges_[i+1] = this->neighbors_.size() + std::distance(first, last);

        // allcoate if the current size is not enough.
        // XXX
        // This invalidates iterators because of the re-allocation.
        // So the `ranges_` cannot have iterators.
        // The elements of `ranges_` must be indices.
        this->neighbors_.reserve(neighbors_.size() + std::distance(first, last));

        std::copy(first, last, std::back_inserter(this->neighbors_));
        return ;
    }

    range_type operator[](const std::size_t i) const noexcept
    {
        return range_type{
            this->neighbors_.begin() + this->ranges_[i],
            this->neighbors_.begin() + this->ranges_[i+1]
        };
    }
    range_type at(const std::size_t i) const
    {
        return range_type{
            this->neighbors_.begin() + this->ranges_.at(i),
            this->neighbors_.begin() + this->ranges_.at(i+1)
        };
    }

  private:
    container_type           neighbors_;
    std::vector<std::size_t> ranges_;
};

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
