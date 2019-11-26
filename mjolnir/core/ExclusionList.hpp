#ifndef MJOLNIR_CORE_EXCLUSION_LIST_HPP
#define MJOLNIR_CORE_EXCLUSION_LIST_HPP
#include <mjolnir/core/IgnoreMolecule.hpp>
#include <mjolnir/core/IgnoreGroup.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/range.hpp>
#include <mjolnir/util/logger.hpp>
#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

namespace mjolnir
{

// In some cases, like excluded volume interaction between particles that
// are connected by a bond, some particular pairs of particles are excluded
// from interacting pairs.
// This class constructs a list that contains a list of pairs that are excluded
// from interacting pairs using information in a topology.
template<typename traitsT>
class ExclusionList
{
  public:
    using traits_type          = traitsT;
    using topology_type        = Topology;
    using molecule_id_type     = topology_type::molecule_id_type;
    using group_id_type        = topology_type::group_id_type;
    using connection_kind_type = topology_type::connection_kind_type;
    using system_type          = System<traits_type>;

    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using ignore_topology_type = // convert from map to vector to make loop fast
        std::vector<std::pair<connection_kind_type, std::size_t>>;

  public:

    ExclusionList(
        const std::map<connection_kind_type, std::size_t>& ignore_top,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : ignore_molecule_(std::move(ignore_mol)),
          ignore_group_   (std::move(ignore_grp)),
          ignore_topology_(ignore_top.begin(), ignore_top.end())
    {}
    ~ExclusionList() = default;
    ExclusionList(const ExclusionList&) = default;
    ExclusionList(ExclusionList&&)      = default;
    ExclusionList& operator=(const ExclusionList&) = default;
    ExclusionList& operator=(ExclusionList&&)      = default;

    // check an interaction exists between i-th and j-th particles.
    //
    // It searches a small list of indices and check j is found in a exclusion-
    // list of i-th particle.
    bool is_excluded(const std::size_t i, const std::size_t j) const
    {
        // assuming all the lists are enough small (< 20 or so)
        {
            const auto ign_grps = this->ignored_grp_of(this->grp_ids_[i]);
            const auto grp_of_j = this->grp_ids_[j];
            if(std::binary_search(ign_grps.begin(), ign_grps.end(), grp_of_j))
            {
                // if found, the pair is ignored. return true.
                return true;
            }
        }

        // check molecule ids...
        const auto mol_of_j = this->mol_ids_[j];
        for(const auto& ignoring_mol : this->ignored_mol_of(this->mol_ids_[i]))
        {
            // already sorted like ignoring_mol = [4,5,6]
            if     (mol_of_j <  ignoring_mol) {break;}
            else if(mol_of_j == ignoring_mol) {return true;}
        }
        // check distance on topology...
        for(const auto& ignoring_idx : this->ignored_idx_of(i))
        {
            if     (ignoring_idx >  j) {break;}
            else if(ignoring_idx == j) {return true;}
        }
        return false;
    }

    void make(const system_type& sys)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        const auto& topol = sys.topology();
        const std::size_t N = sys.size();

        {
            // all groups defined in `ignore.group` table. check mistakes
            auto groups_defined = this->ignore_group_.all_groups();

            // count groups and assign indices
            this->grp_list_.clear();
            for(std::size_t i=0; i<N; ++i)
            {
                if(std::find(this->grp_list_.begin(), this->grp_list_.end(),
                             topol.group_of(i)) == this->grp_list_.end())
                {
                    // if not found, push it to the list
                    this->grp_list_.push_back(topol.group_of(i));
                    MJOLNIR_LOG_INFO("group ", topol.group_of(i), " found. ",
                                     this->grp_list_.size(), "-th group is ",
                                     topol.group_of(i));
                }

                const auto found = std::find(groups_defined.begin(),
                        groups_defined.end(), topol.group_of(i));
                if(found != groups_defined.end())
                {
                    groups_defined.erase(found);
                }
            }
            MJOLNIR_LOG_INFO("all groups are found {", this->grp_list_, "}.");
            if(!groups_defined.empty())
            {
                MJOLNIR_LOG_WARN("unknown group is specified in `ignore.group`");
                for(const auto& g : groups_defined)
                {
                    MJOLNIR_LOG_WARN("- ", g);
                }
            }
        }

        // Here, the list of all groups is constructed.
        // Next, construct grp_ids_. assign corresponding group-ids
        this->grp_ids_.resize(N);
        for(std::size_t i=0; i<N; ++i)
        {
            const auto found = std::find(grp_list_.begin(), grp_list_.end(),
                                         topol.group_of(i));
            assert(found != grp_list_.end()); // should be found.
            this->grp_ids_.at(i) = std::distance(grp_list_.begin(), found);
            MJOLNIR_LOG_DEBUG("particle ", i, " belongs to ", *found, " and "
                              "the corresponding group_idx is ", grp_ids_.at(i),
                              " and the value is ", grp_list_.at(grp_ids_.at(i)));
        }

        // Next, construct exclusion list for group.
        {
            const std::size_t Ngrps = this->grp_list_.size();

            std::size_t idx = 0; // the current index in the exclusion list
            for(std::size_t i=0; i<Ngrps; ++i)
            {
                const std::size_t first = idx; // first index of the sub-list

                const auto& grp_name_i = this->grp_list_.at(i);
                for(std::size_t j=0; j<Ngrps; ++j)
                {
                    const auto& grp_name_j = this->grp_list_.at(j);
                    if(this->is_ignored_group(grp_name_i, grp_name_j))
                    {
                        MJOLNIR_LOG_INFO("group ", grp_name_i, " and ",
                                         grp_name_j, " ignores each other");
                        this->ignored_grps_.push_back(j);
                        ++idx;
                    }
                }
                const auto beg = ignored_grps_.begin();
                // sort this sub-range [first, idx)
                std::sort(beg + first, beg + idx);
                this->grp_ranges_.emplace_back(first, idx);
            }
        }

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
                    if(this->is_ignored_molecule(i, j))
                    {
                        MJOLNIR_LOG_INFO("molecule ", i, " and molecule ", j,
                                         " will ignore each other");
                        this->ignored_mols_.push_back(j);
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
                std::vector<std::size_t> ignored_particles{i}; // ignore itself
                for(const auto& connection : this->ignore_topology_)
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

    // ------------------------------------------------------------------------
    // utilities

    ignore_topology_type const& ignore_topology() const noexcept
    {
        return ignore_topology_;
    }
    bool is_ignored_molecule(
            const molecule_id_type& i, const molecule_id_type& j) const noexcept
    {
        return ignore_molecule_.is_ignored(i, j);
    }
    bool is_ignored_group(
            const group_id_type& i, const group_id_type& j) const noexcept
    {
        return ignore_group_.is_ignored(i, j);
    }

    const char* ignored_molecule_type() const noexcept
    {
        return this->ignore_molecule_.name();
    }
    std::map<group_id_type, std::vector<group_id_type>> const&
    ignored_groups() const noexcept
    {
        return this->ignore_group_.ignores();
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
    range<typename std::vector<std::size_t>::const_iterator>
    ignored_grp_of(const std::size_t i) const noexcept
    {
        return range<typename std::vector<std::size_t>::const_iterator>{
            this->ignored_grps_.begin() + this->grp_ranges_[i].first,
            this->ignored_grps_.begin() + this->grp_ranges_[i].second
        };
    }

  private:

    ignore_molecule_type ignore_molecule_;
    ignore_group_type    ignore_group_;
    ignore_topology_type ignore_topology_;

    // It contains the same infromation as {topol.nodes_.molecule_id};
    std::vector<molecule_id_type> mol_ids_;

    // Here, ExclusionList counts a number of groups and assign a unique index
    // for each group. `grp_ids` contains the index, not std::string itself.
    // It would be changed by calling ExclusionList::make() because it counts
    // groups and and assigns IDs.
    std::vector<group_id_type>    grp_list_; // list of the group names
    std::vector<std::size_t>      grp_ids_;  // the index of the above list

    // ignored mol_id...
    std::vector<std::size_t> ignored_mols_;
    std::vector<std::pair<std::size_t, std::size_t>> mol_ranges_;

    // grp_idx...
    std::vector<std::size_t> ignored_grps_;
    std::vector<std::pair<std::size_t, std::size_t>> grp_ranges_;

    // normal particle idx
    std::vector<std::size_t> ignored_idxs_;
    std::vector<std::pair<std::size_t, std::size_t>> idx_ranges_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ExclusionList<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ExclusionList<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ExclusionList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ExclusionList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_NEIGHBOR_LIST_HPP
