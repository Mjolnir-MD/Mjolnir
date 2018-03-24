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
    typedef topology_type::chain_id_type chain_id_type;
    typedef topology_type::connection_kind_type connection_kind_type;

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
        // assuming both list is enough small (< 20 or so)
        const auto chain_of_j = this->chain_ids_[j];
        for(const auto& ignoring_chn : this->ignored_chn_of(this->chain_ids_[i]))
        {
            // already sorted like ignoring_chn = [4,5,6]
            if     (chain_of_j <  ignoring_chn) {break;}
            else if(chain_of_j == ignoring_chn) {return true;}
        }
        for(const auto& ignoring_idx : this->ignored_idxs_of(i))
        {
            if     (ignoring_idx >  j) {break;}
            else if(ignoring_idx == j) {return true;}
        }
        return false;
    }

    template<typename PotentialT>
    void make(const System& sys, const PotentialT& pot)
    {
        const auto& topol = sys.topology();
        const std::size_t N = sys.size();

        // copy chain_ids from topol to this
        const std::size_t Nchain = topol.number_of_chains();
        this->chain_ids_.resize(N);
        for(std::size_t i=0; i<N; ++i)
        {
            this->chain_ids_[i] = topol.chain_of(i);
        }

        // make ignored_chain_idxs
        {
            std::size_t idx = 0;
            for(std::size_t i=0; i<Nchain; ++i)
            {
                const std::size_t first = idx;
                for(std::size_t j=0; j<Nchain; ++j)
                {
                    if(pot.is_ignored_chain(i, j))
                    {
                        this->ignored_chns_.push_back(i);
                        ++idx;
                    }
                }
                const auto beg = ignored_chns_.begin();
                std::sort(beg + first, beg + idx);
                this->chn_ranges_.emplace_back(first, idx);
            }
        }

        // make ignored_particle_idxs
        // excluded_connection := pair{connection kind, distance}
        {
            std::size_t idx = 0;
            for(std::size_t i=0; i<N; ++i)
            {
                const std::size_t first = idx;
                {
                    const std::size_t dist = pot.ignored_bonds();
                    for(const auto j : topol.list_adjacent_within(
                                i, dist, connection_kind_type::bond))
                    {
                        this->ignored_idxs_.push_back(j);
                        ++idx;
                    }
                }
                {
                    const std::size_t dist = pot.ignored_contacts();
                    for(const auto j : topol.list_adjacent_within(
                                i, dist, connection_kind_type::contact))
                    {
                        this->ignored_idxs_.push_back(j);
                        ++idx;
                    }
                }
                const auto beg = ignored_idxs_.begin();
                std::sort(beg + first, beg + idx);
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
    ignored_chn_of(const std::size_t i) const noexcept
    {
        return range_type{
            this->ignored_chns_.begin() + this->chn_ranges_[i].first,
            this->ignored_chns_.begin() + this->chn_ranges_[i].second
        };
    }

  private:

    std::vector<chain_id_type> chain_ids_; // same as {topol.nodes_.chain_id};

    // ignored chains...
    std::vector<std::size_t> ignored_chns_;
    std::vector<std::pair<std::size_t, std::size_t>> chn_ranges_;

    std::vector<std::size_t> ignored_idxs_;
    std::vector<std::pair<std::size_t, std::size_t>> idx_ranges_;
};

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
