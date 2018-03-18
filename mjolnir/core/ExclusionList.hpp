#ifndef MJOLNIR_EXCLUSION_LIST_HPP
#define MJOLNIR_EXCLUSION_LIST_HPP
#include <mjolnir/util/range.hpp>
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

namespace mjolnir
{

class ExclusionList
{
  public:
    typedef StructureTopology topology_type;
    typedef topology_type::group_id_type        group_id_type;
    typedef topology_type::connection_name_type connection_name_type;
    typedef std::pair<connection_name_type, std::size_t> ignored_connection_type;

  public:

    ExclusionList()
      : ignored_group([](group_id_type, group_id_type) -> bool {return false;})
    {}
    ~ExclusionList() = default;
    ExclusionList(const ExclusionList&) = default;
    ExclusionList(ExclusionList&&)      = default;
    ExclusionList& operator=(const ExclusionList&) = default;
    ExclusionList& operator=(ExclusionList&&)      = default;

    template<typename BinaryPredicate>
    void set_group_exclusion(BinaryPredicate&& pred)
    {
        // if group i and j should be ignored, return true.
        this->ignored_group = std::forward<BinaryPredicate>(pred);
        return;
    }

    /*! @brief predicates whether i and j should be ignored
     *  @param i particle id
     *  @param j particle id */
    bool is_excluded(const std::size_t i, const std::size_t j) const
    {
        // XXX std::function::operator() checks whether functor null or not
        //     everytime called. consider making list?
        if(ignored_group(this->group_ids_[i], this->group_ids_[j]))
        {
            return true;
        }
        for(const auto ignoring_idx : this->ignored_idxs_of(i))
        {
            if(ignoring_idx == j) {return true;}
        }
        return false;
    }

    void make(const StructureTopology& topol,
              const std::vector<ignored_connection_type>& exclude_connection)
    {
        const std::size_t N = topol.size();

        // copy group ids from topol to this
        this->group_ids_.resize(N);
        for(std::size_t i=0; i<N; ++i)
        {
            this->group_ids_[i] = topol.group_of(i);
        }

        // make ignored_particle_idxs
        std::size_t idx = 0;
        for(std::size_t i=0; i<N; ++i)
        {
            const std::size_t first = idx;
            for(const auto& connection_kind : exclude_connection)
            {
                const connection_name_type& name = connection_kind.first;
                const std::size_t           dist = connection_kind.second;

                for(const auto cnct : topol.list_adjacent_within(i, dist, name))
                {
                    this->ignored_idxs_.push_back(cnct);
                    ++idx;
                }
            }
            this->idx_ranges_.emplace_back(first, idx);
        }
        return;
    }

  private:

    range<typename std::vector<std::size_t>::const_iterator>
    ignored_idx_of(const std::size_t i) const noexcept
    {
        return range_type{
            this->ignored_idxs_.begin() + this->idx_ranges_[i].first,
            this->ignored_idxs_.begin() + this->idx_ranges_[i].second
        };
    }

  private:

    std::vector<group_id_type> group_ids_; // same as {topol.nodes_.group_id};
    std::function<bool(group_id_type, group_id_type)> ignored_group;

    std::vector<std::size_t> ignored_idxs_;
    std::vector<std::pair<std::size_t, std::size_t>> idx_ranges_;
};

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
