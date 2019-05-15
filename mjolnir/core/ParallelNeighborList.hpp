#ifndef MJOLNIR_CORE_PARALLEL_NEIGHBOR_LIST_HPP
#define MJOLNIR_CORE_PARALLEL_NEIGHBOR_LIST_HPP
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// This class is introduced to enable parallel construction of a neighbor list.
// ParallelNeighborList::add_list_for is thread-safe.
//
// The difference from Normal NeighborList is that it first receives the
// maximum number of interacting partners. For example, let's consider
// particle 1 has 5 partners, p2 has 6, p3 has 4, then the number is 6.
// In that case, to construct neighbor list in parallel, thus we need at least 6
// elements dedicated for each particle. Also, to avoid false-sharing, it is
// good to align the sub-range of partners by cache-alignment.
//
// So `ParallelNeighborList::resize` *must* be called before calling
// `add_list_for`. Without calling that, it does not have enough memory
// to access to the sub-range.
template<typename parameterT>
class ParallelNeighborList
{
  public:
    using parameter_type = parameterT;
    using neighbor_type  = neighbor_element<parameter_type>;
    using container_type = std::vector<neighbor_type>;
    using range_type     = range<typename container_type::const_iterator>;

  public:

    ParallelNeighborList()  = default;
    ~ParallelNeighborList() = default;
    ParallelNeighborList(const ParallelNeighborList&) = default;
    ParallelNeighborList(ParallelNeighborList&&)      = default;
    ParallelNeighborList& operator=(const ParallelNeighborList&) = default;
    ParallelNeighborList& operator=(ParallelNeighborList&&)      = default;

    void clear()
    {
        this->neighbors_.clear();
        this->sizes_    .clear();
        return;
    }

    void resize(const std::size_t Nparticle, const std::size_t max_neighbor)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("Nparticle                = ", Nparticle);
        MJOLNIR_LOG_INFO("max_neighbor             = ", max_neighbor);

        constexpr std::size_t cache_line = 64; // bytes
        const auto max_partner_size = max_neighbor * sizeof(neighbor_type);

        MJOLNIR_LOG_INFO("max_partner_size         = ", max_partner_size);

        if(max_partner_size % cache_line == 0)
        {
            this->max_partner_per_particle_ = max_neighbor;
        }
        else
        {
            this->max_partner_per_particle_ = (cache_line *
                (max_partner_size / cache_line + 1)) / sizeof(neighbor_type);
        }

        MJOLNIR_LOG_INFO("max_partner_per_particle = ", max_partner_per_particle_);

        this->neighbors_.resize(Nparticle * max_partner_per_particle_);
        this->sizes_    .resize(Nparticle, /* default value is */ 0);
        return;
    }

    // assign the partners in range [first, last) as the partner of particle i.
    template<typename Iterator>
    void add_list_for(const std::size_t i, Iterator first, Iterator last)
    {
        static_assert(std::is_same<
            typename std::iterator_traits<Iterator>::value_type, neighbor_type
            >::value, "The iterator must points neighbor type.");

        assert(i < sizes_.size());
        assert(std::size_t(std::distance(first, last)) < max_partner_per_particle_);

        assert(this->neighbors_.size() >= max_partner_per_particle_ * (i+1));
        const std::size_t offset = this->max_partner_per_particle_ * i;

        this->sizes_[i] = std::distance(first, last);
        std::copy(first, last, this->neighbors_.begin() + offset);
        return ;
    }

    range_type operator[](const std::size_t i) const noexcept
    {
        const std::size_t offset = this->max_partner_per_particle_ * i;
        return range_type{
            this->neighbors_.begin() + offset,
            this->neighbors_.begin() + offset + this->sizes_[i]
        };
    }
    range_type at(const std::size_t i) const
    {
        const std::size_t offset = this->max_partner_per_particle_ * i;
        return range_type{
            this->neighbors_.begin() + offset,
            this->neighbors_.begin() + offset + this->sizes_.at(i)
        };
    }

  private:
    std::size_t              max_partner_per_particle_;
    container_type           neighbors_;
    std::vector<std::size_t> sizes_;
};

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
