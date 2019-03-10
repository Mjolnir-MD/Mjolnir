#ifndef MJOLNIR_EXCLUSION_LIST_HPP
#define MJOLNIR_EXCLUSION_LIST_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/range.hpp>
#include <mjolnir/util/logger.hpp>
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

namespace mjolnir
{

class ExclusionList
{
  public:
    using topology_type        = Topology;
    using molecule_id_type     = topology_type::molecule_id_type;
    using connection_kind_type = topology_type::connection_kind_type;

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
        const auto mol_of_j = this->mol_ids_[j];
        for(const auto& ignoring_mol : this->ignored_mol_of(this->mol_ids_[i]))
        {
            // already sorted like ignoring_mol = [4,5,6]
            if     (mol_of_j <  ignoring_mol) {break;}
            else if(mol_of_j == ignoring_mol) {return true;}
        }
        for(const auto& ignoring_idx : this->ignored_idx_of(i))
        {
            if     (ignoring_idx >  j) {break;}
            else if(ignoring_idx == j) {return true;}
        }
        return false;
    }

    template<typename traitsT, typename PotentialT>
    void make(const System<traitsT>& sys, const PotentialT& pot)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential = ", pot.name());

        const auto& topol = sys.topology();
        const std::size_t N = sys.size();

        // copy molecule_ids from topol to this
        this->mol_ids_.resize(N);
        for(std::size_t i=0; i<N; ++i)
        {
            this->mol_ids_[i] = topol.molecule_of(i);
            MJOLNIR_LOG_DEBUG("particle ", i, " is belonging molecule ",
                              topol.molecule_of(i));
        }

        // make ignored_molecule_idxs
        {
            const std::size_t Nmols = topol.number_of_molecules();
            std::size_t idx = 0;
            for(std::size_t i=0; i<Nmols; ++i)
            {
                const std::size_t first = idx;
                for(std::size_t j=0; j<Nmols; ++j)
                {
                    if(pot.is_ignored_molecule(i, j))
                    {
                        MJOLNIR_LOG_INFO("molecule ", i, " and molecule ", j,
                            " will ignore each other on ", pot.name());
                        this->ignored_mols_.push_back(i);
                        ++idx;
                    }
                }
                const auto beg = ignored_mols_.begin();
                std::sort(beg + first, beg + idx);
                this->mol_ranges_.emplace_back(first, idx);
            }
        }

        // make ignored_particle_idxs
        // excluded_connection := pair{connection kind, distance}
        {
            std::size_t idx = 0;
            for(std::size_t i=0; i<N; ++i)
            {
                const std::size_t first = idx;
                std::vector<std::size_t> ignored_particles;
                for(const auto& connection : pot.ignore_within())
                {
                    const std::size_t dist = connection.second;
                    for(const auto j :
                        topol.list_adjacent_within(i, dist, connection.first))
                    {
                        ignored_particles.push_back(j);
                    }
                }
                std::sort(ignored_particles.begin(), ignored_particles.end());
                const auto last = std::unique(ignored_particles.begin(),
                                              ignored_particles.end());
                ignored_particles.erase(last, ignored_particles.end());
                MJOLNIR_LOG_INFO("particle ", i, " ignores ", ignored_particles);

                for(const auto j : ignored_particles)
                {
                    this->ignored_idxs_.push_back(j);
                    ++idx;
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
        return range<typename std::vector<std::size_t>::const_iterator>{
            this->ignored_idxs_.begin() + this->idx_ranges_[i].first,
            this->ignored_idxs_.begin() + this->idx_ranges_[i].second
        };
    }
    range<typename std::vector<std::size_t>::const_iterator>
    ignored_mol_of(const std::size_t i) const noexcept
    {
        return range<typename std::vector<std::size_t>::const_iterator>{
            this->ignored_mols_.begin() + this->mol_ranges_[i].first,
            this->ignored_mols_.begin() + this->mol_ranges_[i].second
        };
    }

  private:

    // It contains the same infromation as {topol.nodes_.molecule_id};
    std::vector<molecule_id_type> mol_ids_;

    // ignored mol_id...
    std::vector<std::size_t> ignored_mols_;
    std::vector<std::pair<std::size_t, std::size_t>> mol_ranges_;

    std::vector<std::size_t> ignored_idxs_;
    std::vector<std::pair<std::size_t, std::size_t>> idx_ranges_;
};

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
