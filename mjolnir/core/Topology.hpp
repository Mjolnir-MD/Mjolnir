#ifndef MJOLNIR_CORE_TOPOLOGY_HPP
#define MJOLNIR_CORE_TOPOLOGY_HPP
#include <mjolnir/util/throw_exception.hpp>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace mjolnir
{

// XXX Topology class manages the connections between beads.
//
// NOTE: connection kind "bond" is used to determine the molecules.
//       Other kind of connections also can be defined (and overlayed
//       if nessesary). In that case, you can traverse the graph by using
//       your connection kind also, but it will not be considered when
//       it detemines "molecule"s.
//
// Most of the time, "molecule" defined by this class is equivalent to "chain".
// But calling "molecule" as "chain" is misreading because there can be a
// circular molecule or tiny molecule that only has a few number of beads.
// For example, lipid membrane is composed of a large number of lipid molecules
// and people useally does not call each lipid as "chain".
//
// Generally, traversing graph takes time rather than finding value from
// contiguous container. So there are ExclusionList class that stores the
// information for NeighborLists (usually, nonlocal potentials are ignored if
// the two particles are connected by covarent bond). While the simulation,
// ExclusionLists would be called many time, but Topology would only be called
// to construct ExclusionList at the beginning ofjthe simulation.
//
// Since the members of this class are rarely called, in most cases, the speed
// of these methods does not affect to the overall efficiency.class Topology
class Topology
{
  public:

    using molecule_id_type     = std::size_t;
    using group_id_type        = std::string;
    using name_type            = std::string;
    using connection_kind_type = std::string;
    using edge_type            = std::pair<std::size_t, connection_kind_type>;

    static constexpr molecule_id_type uninitialized() noexcept
    {
        return std::numeric_limits<molecule_id_type>::max();
    }

    struct node
    {
        std::size_t molecule_id;
        std::vector<edge_type> adjacents;
    };

  public:

    Topology()  = default;
    ~Topology() = default;
    Topology(const Topology&) = default;
    Topology(Topology&&)      = default;
    Topology& operator=(const Topology&) = default;
    Topology& operator=(Topology&&)      = default;

    explicit Topology(const std::size_t N)
        : num_molecules_(1), nodes_(N, node{uninitialized(), {}}),
          names_(N, "uninitialized"), groups_(N, "uninitialized")
    {}

    void        clear()                {return nodes_.clear();}
    bool        empty() const noexcept {return nodes_.empty();}
    std::size_t size()  const noexcept {return nodes_.size();}

    name_type&       name_of(const std::size_t i)       {return this->names_.at(i);}
    name_type const& name_of(const std::size_t i) const {return this->names_.at(i);}
    name_type&       name_of(const std::size_t i, const std::nothrow_t&)       noexcept {return this->names_[i];}
    name_type const& name_of(const std::size_t i, const std::nothrow_t&) const noexcept {return this->names_[i];}

    molecule_id_type  molecule_of(const std::size_t i) const
    {return nodes_.at(i).molecule_id;}
    molecule_id_type& molecule_of(const std::size_t i)
    {return nodes_.at(i).molecule_id;}

    molecule_id_type  molecule_of(const std::size_t i, const std::nothrow_t&) const
    {return nodes_[i].molecule_id;}
    molecule_id_type& molecule_of(const std::size_t i, const std::nothrow_t&)
    {return nodes_[i].molecule_id;}

    group_id_type&       group_of(const std::size_t i)       {return this->groups_.at(i);}
    group_id_type const& group_of(const std::size_t i) const {return this->groups_.at(i);}
    group_id_type&       group_of(const std::size_t i, const std::nothrow_t&)       noexcept {return this->groups_[i];}
    group_id_type const& group_of(const std::size_t i, const std::nothrow_t&) const noexcept {return this->groups_[i];}

    void add_connection  (const std::size_t i, const std::size_t j,
                          const connection_kind_type& kind);
    void erase_connection(const std::size_t i, const std::size_t j,
                          const connection_kind_type& kind);
    bool has_connection  (const std::size_t i, const std::size_t j,
                          const connection_kind_type& kind);
    // has_connection returns true if i == j.

    std::vector<std::size_t>
    list_adjacent_within(const std::size_t node_idx, const std::size_t dist,
                         const connection_kind_type& kind) const;

    std::vector<connection_kind_type>
    list_connections_between(const std::size_t i, const std::size_t j) const;

    //! reset molecule_id of all the particles
    void construct_molecules();
    void resize(const std::size_t N) {nodes_.resize(N); return;}

    std::size_t number_of_molecules() const noexcept {return this->num_molecules_;}


  private:

    void list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist,
        const connection_kind_type& kind, std::vector<std::size_t>& out) const
    {
        if(dist == 0)
        {
            out.push_back(node_idx);
            return;
        }

        for(auto edge : this->nodes_.at(node_idx).adjacents)
        {
            if(edge.second != kind)
            {
                continue;
            }
            out.push_back(edge.first);
            this->list_adjacent_within(edge.first, dist-1, kind, out);
        }
        return;
    }

  private:
    // each node corresponds to the particle having the same idx in a system.
    std::size_t                num_molecules_;
    std::vector<node>          nodes_;
    std::vector<name_type>     names_;
    std::vector<group_id_type> groups_;
};

inline void Topology::add_connection(
        const std::size_t i, const std::size_t j,
        const connection_kind_type& kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw_exception<std::out_of_range>(
            "mjolnir::Topology::add_connection: size of nodes = ",
            nodes_.size(), ", i = ", i, ", j = ", j);
    }

    {
        const edge_type edge = std::make_pair(j, kind);

        auto& adjs = nodes_[i].adjacents;
        if(std::find(adjs.cbegin(), adjs.cend(), edge) == adjs.cend())
        {
            adjs.push_back(edge);
        }
    }
    {
        const edge_type edge = std::make_pair(i, kind);

        auto& adjs = nodes_[j].adjacents;
        if(std::find(adjs.cbegin(), adjs.cend(), edge) == adjs.cend())
        {
            adjs.push_back(edge);
        }
    }
    return ;
}

inline void Topology::erase_connection(
        const std::size_t i, const std::size_t j,
        const connection_kind_type& kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw_exception<std::out_of_range>(
            "mjolnir::Topology::erase_connection: size of nodes = ",
            nodes_.size(), ", i = ", i, ", j = ", j);
    }

    {
        auto& adjs = nodes_[i].adjacents;
        const auto found = std::find(
                adjs.begin(), adjs.end(), std::make_pair(j, kind));
        if(found != adjs.end())
        {
            adjs.erase(found);
        }
    }
    {
        auto& adjs = nodes_[j].adjacents;
        const auto found = std::find(
                adjs.begin(), adjs.end(), std::make_pair(i, kind));
        if(found != adjs.end())
        {
            adjs.erase(found);
        }
    }
    return ;
}

inline bool Topology::has_connection(
        const std::size_t i, const std::size_t j,
        const connection_kind_type& kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw_exception<std::out_of_range>(
            "mjolnir::Topology::has_connection: size of nodes = ",
            nodes_.size(), ", i = ", i, ", j = ", j);
    }
    if(i == j) {return true;} // XXX

    const auto& adjs = nodes_[i].adjacents;
    const auto found =
        std::find(adjs.cbegin(), adjs.cend(), std::make_pair(j, kind));
    return found != adjs.cend();
}

inline std::vector<std::size_t>
Topology::list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist,
        const connection_kind_type& kind) const
{
    std::vector<std::size_t> retval = {node_idx};
    if(dist == 0) {return retval;}

    for(auto edge : this->nodes_.at(node_idx).adjacents)
    {
        if(edge.second != kind)
        {
            continue;
        }
        retval.push_back(edge.first);
        this->list_adjacent_within(edge.first, dist - 1, kind, retval);
    }
    std::sort(retval.begin(), retval.end());
    const auto new_end = std::unique(retval.begin(), retval.end());
    retval.erase(new_end, retval.end());
    return retval;
}

inline std::vector<typename Topology::connection_kind_type>
Topology::list_connections_between(const std::size_t i, const std::size_t j) const
{
    std::vector<connection_kind_type> connections;
    for(const auto& edge : this->nodes_.at(i).adjacents)
    {
        if(edge.first == j)
        {
            connections.push_back(edge.second);
        }
    }
    return connections;
}

inline void
Topology::construct_molecules()
{
    if(this->nodes_.empty()){return;}
    // reset molecule-ids
    for(auto& node : nodes_)
    {
        node.molecule_id = uninitialized();
    }

    molecule_id_type next_molecule_id = 0;
    for(auto& node : nodes_)
    {
        // look all the adjacents and update node by their molecule id
        for(const auto& edge : node.adjacents)
        {
            // ignore all the edges that are not bonds
            if(edge.second != "bond") {continue;}

            // if an adjacent that is connected by a "bond",
            // it should be in the same molecule.
            const auto& adj = this->nodes_.at(edge.first);
            if(adj.molecule_id != uninitialized())
            {
                node.molecule_id = adj.molecule_id;
                break;
            }
        }
        // if the molecule_id of this node is still not initialized, that means
        // this node belongs to a new molecule.
        if(node.molecule_id == uninitialized())
        {
            node.molecule_id = next_molecule_id++;
        }
    }
    this->num_molecules_ = next_molecule_id;
    return;
}

} // mjolnir
#endif// MJOLNIR_STRUCTURE_TOPOLOGY_H
