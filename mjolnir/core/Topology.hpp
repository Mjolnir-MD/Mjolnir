#ifndef MJOLNIR_CORE_TOPOLOGY_H
#define MJOLNIR_CORE_TOPOLOGY_H
#include <mjolnir/util/throw_exception.hpp>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace mjolnir
{

// XXX Topology class manages the connections between beads.
// Generally, traversing graph takes time rather than finding value from
// contiguous container. So there are ExclusionList class that stores the
// connection information for NeighborLists (usually, nonlocal potentials are
// ignored if the two particles are connected by covarent bond).
// Since the members of this class are rarely called, in most cases, the speed
// of these methods does not affect to the overall efficiency.
class Topology
{
  public:

    using chain_id_type        = std::size_t;
    using connection_kind_type = std::string;
    using edge_type = std::pair<std::size_t, connection_kind_type>;

    // there is one `reserved` connection kind that has special meaning.
    // "bond" is used to determine chains.

    static constexpr chain_id_type uninitialized =
        std::numeric_limits<chain_id_type>::max();

    // each node corresponds to the particle having the same idx in a system.
    struct node
    {
        std::size_t chain_id;
        std::string identifier;
        std::vector<edge_type> adjacents;
    };

  public:

    Topology()  = default;
    ~Topology() = default;
    Topology(const Topology&) = default;
    Topology(Topology&&)      = default;
    Topology& operator=(const Topology&) = default;
    Topology& operator=(Topology&&)      = default;

    Topology(const std::size_t N)
        : nodes_(N, {uninitialized, "uninitialized", {}})
    {}

    void        clear()                {return nodes_.clear();}
    bool        empty() const noexcept {return nodes_.empty();}
    std::size_t size()  const noexcept {return nodes_.size();}

    std::string&       identifier_of(const std::size_t i)
    {return this->nodes_.at(i).identifier;}
    std::string const& identifier_of(const std::size_t i) const
    {return this->nodes_.at(i).identifier;}

    chain_id_type  chain_of(const std::size_t idx) const
    {return nodes_.at(idx).chain_id;}
    chain_id_type& chain_of(const std::size_t idx)
    {return nodes_.at(idx).chain_id;}

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
    list_connections_between(const std::size_t i, const std::size_t j,
                             const connection_kind_type& kind) const;

    //! reset chain_id of all the particles
    void construct_chains();
    void resize(const std::size_t N) {nodes_.resize(N); return;}

    std::size_t number_of_chains() const noexcept {return this->num_chains_;}


  private:

    void list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist,
        const connection_kind_type& kind, std::vector<std::size_t>& out) const
    {
        if(dist == 0)
        {
            if(std::find(out.cbegin(), out.cend(), node_idx) != out.cend())
            {
                out.push_back(node_idx);
            }
            return;
        }

        for(auto edge : this->nodes_.at(node_idx).adjacents)
        {
            if(edge.second != kind ||
               std::find(out.cbegin(), out.cend(), edge.first) != out.cend())
            {
                continue;
            }
            out.push_back(edge.first);
            this->list_adjacent_within(edge.first, dist-1, kind, out);
        }
        return;
    }

  private:
    std::size_t  num_chains_;
    std::vector<node> nodes_;
};
constexpr typename Topology::chain_id_type Topology::uninitialized;

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
        if(edge.second != kind || std::find(
           retval.cbegin(), retval.cend(), edge.first) != retval.cend())
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
Topology::list_connections_between(const std::size_t i, const std::size_t j,
                                   const connection_kind_type& kind) const
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
Topology::construct_chains()
{
    if(this->nodes_.empty()){return;}
    for(auto& node : nodes_)
    {
        node.chain_id = uninitialized;
    }

    chain_id_type next_chain_id = 0;
    for(auto& node : nodes_)
    {
        for(const auto& edge : node.adjacents)
        {
            if(edge.second != "bond") {continue;}
            const auto& adj = this->nodes_.at(edge.first);
            if(adj.chain_id != uninitialized)
            {
                node.chain_id = adj.chain_id;
                break;
            }
        }
        if(node.chain_id == uninitialized)
        {
            node.chain_id = next_chain_id++;
        }
    }
    this->num_chains_ = ++next_chain_id;
    return;
}

} // mjolnir
#endif// MJOLNIR_STRUCTURE_TOPOLOGY_H