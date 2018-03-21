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
    typedef topology_type::group_id_type group_id_type;

  public:

    ExclusionList() = default;
    ~ExclusionList() = default;
    ExclusionList(const ExclusionList&) = default;
    ExclusionList(ExclusionList&&)      = default;
    ExclusionList& operator=(const ExclusionList&) = default;
    ExclusionList& operator=(ExclusionList&&)      = default;

    /*! @brief predicates whether i and j should be ignored
     *  @param i particle id
     *  @param j particle id */
    bool is_excluded(const std::size_t i, const std::size_t j) const
    {
        const group_id_type i_grp = this->group_ids_[i];
        const group_id_type j_grp = this->group_ids_[j];

        for(const auto ignoring_grp : this->ignored_grp_of(i_grp))
        {
            if(ignoring_grp == j_grp) {return true;}
        }
        for(const auto ignoring_idx : this->ignored_idxs_of(i))
        {
            if(ignoring_idx == j) {return true;}
        }
        return false;
    }

    template<typename PotentialT>
    void make(const System& sys, const PotentialT& pot)
    {
        const auto& topol = sys.topology();
        const std::size_t N = sys.size();

        // copy group ids from topol to this
        std::vector<group_id_type> group_kind;
        this->group_ids_.resize(N);
        for(std::size_t i=0; i<N; ++i)
        {
            const auto grp = topol.group_of(i);
            this->group_ids_[i] = grp;

            if(std::find(group_kind.begin(), group_kind.end(), grp) ==
                    group_kind.end())
            {
                group_kind.push_back(grp);
            }
        }

        // make ignored_group_idxs
        {
            std::size_t idx = 0;
            for(std::size_t i=0; i<group_kind.size(); ++i)
            {
                const std::size_t first = idx;
                for(std::size_t j=0; j<group_kind.size(); ++j)
                {
                    if(pot.is_ignored_group(i, j))
                    {
                        this->ignored_grps_.push_back(i);
                        ++idx;
                    }
                }
                this->grp_ranges_.emplace_back(first, idx);
            }
        }

        // make ignored_particle_idxs
        // excluded_connection := pair{connection kind, distance}
        {
            const auto& excluded_connection = pot.excluded_connections();
            std::size_t idx = 0;
            for(std::size_t i=0; i<N; ++i)
            {
                const std::size_t first = idx;
                for(const auto& connection_kind : excluded_connection)
                {
                    const auto&       name = connection_kind.first;
                    const std::size_t dist = connection_kind.second;

                    for(const auto j :
                            topol.list_adjacent_within(i, dist, name))
                    {
                        this->ignored_idxs_.push_back(j);
                        ++idx;
                    }
                }
                this->idx_ranges_.emplace_back(first, idx);
            }
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
    range<typename std::vector<std::size_t>::const_iterator>
    ignored_grp_of(const std::size_t i) const noexcept
    {
        return range_type{
            this->ignored_grps_.begin() + this->grp_ranges_[i].first,
            this->ignored_grps_.begin() + this->grp_ranges_[i].second
        };
    }

  private:

    std::vector<group_id_type> group_ids_; // same as {topol.nodes_.group_id};

    std::vector<std::size_t> ignored_grps_;
    std::vector<std::pair<std::size_t, std::size_t>> grp_ranges_;

    std::vector<std::size_t> ignored_idxs_;
    std::vector<std::pair<std::size_t, std::size_t>> idx_ranges_;
};

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
