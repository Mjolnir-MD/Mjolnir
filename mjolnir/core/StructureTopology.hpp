#ifndef MJOLNIR_STRUCTURE_TOPOLOGY_H
#define MJOLNIR_STRUCTURE_TOPOLOGY_H
#include <mjolnir/util/throw_exception.hpp>
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace mjolnir
{

class StructureTopology
{
  public:

    typedef std::size_t chain_id_type;
    enum class connection_kind_type
    {
        bond,    //! define chain
        contact, //! does not have effect on the chain (nonbonded contact)
    };

    static constexpr chain_id_type uninitialized() noexcept
    {
        return std::numeric_limits<chain_id_type>::max();
    }

    // each node corresponds to the particle having same idx in a system.
    struct node
    {
        node(): chain_id(uninitialized()) {}
        ~node() = default;
        node(node const&) = default;
        node(node&&)      = default;
        node& operator=(node const&) = default;
        node& operator=(node&&)      = default;

        std::size_t chain_id;
        std::vector<std::pair<std::size_t, connection_kind_type>> adjacents;
    };

  public:

    StructureTopology()  = default;
    ~StructureTopology() = default;
    StructureTopology(const StructureTopology&) = default;
    StructureTopology(StructureTopology&&)      = default;
    StructureTopology& operator=(const StructureTopology&) = default;
    StructureTopology& operator=(StructureTopology&&)      = default;

    StructureTopology(const std::size_t N): nodes_(N){};

    bool        empty() const noexcept {return nodes_.empty();}
    std::size_t size()  const noexcept {return nodes_.size();}

    void resize(const std::size_t N) {nodes_.resize(N); return;}

    void add_connection(const std::size_t i, const std::size_t j,
                        const connection_kind_type kind);

    void erase_connection(const std::size_t i, const std::size_t j);
    void erase_connection(const std::size_t i, const std::size_t j,
                          const connection_kind_type kind);

    //! reset chain_id of all the particles
    void construct_chains();

    bool has_connection(const std::size_t i, const std::size_t j);
    bool has_connection(const std::size_t i, const std::size_t j,
                        const connection_kind_type kind);

    std::vector<std::size_t>
    list_adjacent_within(const std::size_t node_idx, const std::size_t dist
                         ) const;
    std::vector<std::size_t>
    list_adjacent_within(const std::size_t node_idx, const std::size_t dist,
                         const connection_kind_type kind) const;

    chain_id_type number_of_chains() const noexcept {return this->num_chains_;}

    chain_id_type  chain_of(const std::size_t idx) const
    {return nodes_.at(idx).chain_id;}
    chain_id_type& chain_of(const std::size_t idx)
    {return nodes_.at(idx).chain_id;}

  private:
    void
    list_adjacent_within(const std::size_t node_idx, const std::size_t dist,
                         std::vector<std::size_t>& out) const;
    void
    list_adjacent_within(const std::size_t node_idx, const std::size_t dist,
                         const connection_kind_type kind, std::vector<std::size_t>& out
                         ) const;
  private:
    std::size_t  num_chains_;
    std::vector<node> nodes_;
};

inline void StructureTopology::add_connection(
        const std::size_t i, const std::size_t j, connection_kind_type kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw_exception<std::out_of_range>(
            "mjolnir::StructureTopology::add_connection: size of nodes = ",
            nodes_.size(), ", i = ", i, ", j = ", j);
    }

    {
        const auto elem = std::make_pair(j, kind);
        auto& adjs = nodes_[i].adjacents;
        if(std::find(adjs.cbegin(), adjs.cend(), elem) == adjs.cend())
        {
            adjs.push_back(elem);
        }
    }
    {
        const auto elem = std::make_pair(i, kind);
        auto& adjs = nodes_[j].adjacents;
        if(std::find(adjs.cbegin(), adjs.cend(), elem) == adjs.cend())
        {
            adjs.push_back(elem);
        }
    }
    return ;
}

inline void StructureTopology::erase_connection(
        const std::size_t i, const std::size_t j)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw_exception<std::out_of_range>(
            "mjolnir::StructureTopology::erase_connection: size of nodes = ",
            nodes_.size(), ", i = ", i, ", j = ", j);
    }

    {
        auto& adjs = nodes_[i].adjacents;
        const auto found = std::find_if(adjs.begin(), adjs.end(),
            [j](const std::pair<std::size_t, connection_kind_type>& x){
                return x.first == j;
            });
        if(found != adjs.end())
        {
            adjs.erase(found);
        }
    }
    {
        auto& adjs = nodes_[j].adjacents;
        const auto found = std::find_if(adjs.begin(), adjs.end(),
            [i](const std::pair<std::size_t, connection_kind_type>& x){
                return x.first == i;
            });
        if(found != adjs.end())
        {
            adjs.erase(found);
        }
    }
    return ;
}


inline void StructureTopology::erase_connection(
        const std::size_t i, const std::size_t j,
        const connection_kind_type kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw_exception<std::out_of_range>(
            "mjolnir::StructureTopology::erase_connection: size of nodes = ",
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

inline bool StructureTopology::has_connection(
        const std::size_t i, const std::size_t j)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw_exception<std::out_of_range>(
            "mjolnir::StructureTopology::has_connection: size of nodes = ",
            nodes_.size(), ", i = ", i, ", j = ", j);
    }
    if(i == j) {return true;} // XXX

    const auto& adjs = nodes_[i].adjacents;
    const auto found = std::find_if(adjs.cbegin(), adjs.cend(),
        [j](const std::pair<std::size_t, connection_kind_type>& x){
            return x.first == j;
        });

    return found != adjs.cend();
}

inline bool StructureTopology::has_connection(
        const std::size_t i, const std::size_t j,
        const connection_kind_type kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw_exception<std::out_of_range>(
            "mjolnir::StructureTopology::has_connection: size of nodes = ",
            nodes_.size(), ", i = ", i, ", j = ", j);
    }
    if(i == j) {return true;} // XXX

    const auto& adjs = nodes_[i].adjacents;
    const auto found = std::find(adjs.cbegin(), adjs.cend(),
                                 std::make_pair(j, kind));
    return found != adjs.cend();
}

inline std::vector<std::size_t>
StructureTopology::list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist) const
{
    std::vector<std::size_t> retval = {node_idx};
    if(dist == 0) {return retval;}

    for(auto ik : this->nodes_.at(node_idx).adjacents)
    {
        if(std::find(retval.cbegin(), retval.cend(), ik.first) != retval.cend())
        {
            continue; // already assigned. ignore it.
        }
        retval.push_back(ik.first);
        this->list_adjacent_within(ik.first, dist - 1, retval);
    }
    std::sort(retval.begin(), retval.end());
    const auto new_end = std::unique(retval.begin(), retval.end());
    retval.erase(new_end, retval.end());
    return retval;
}

inline std::vector<std::size_t>
StructureTopology::list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist,
        const connection_kind_type kind) const
{
    std::vector<std::size_t> retval = {node_idx};
    if(dist == 0) {return retval;}

    for(auto ik : this->nodes_.at(node_idx).adjacents)
    {
        if(ik.second != kind ||
           std::find(retval.cbegin(), retval.cend(), ik.first) != retval.cend())
        {
            continue;
        }
        retval.push_back(ik.first);
        this->list_adjacent_within(ik.first, dist - 1, kind, retval);
    }
    std::sort(retval.begin(), retval.end());
    const auto new_end = std::unique(retval.begin(), retval.end());
    retval.erase(new_end, retval.end());
    return retval;
}

inline void StructureTopology::list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist,
        std::vector<std::size_t>& out) const
{
    if(dist == 0)
    {
        if(std::find(out.cbegin(), out.cend(), node_idx) != out.cend())
        {
            out.push_back(node_idx);
        }
        return ;
    }

    for(auto ik : this->nodes_.at(node_idx).adjacents)
    {
        if(std::find(out.cbegin(), out.cend(), ik.first) != out.cend())
        {
            continue;
        }
        out.push_back(ik.first);
        this->list_adjacent_within(ik.first, dist - 1, out);
    }
    return;
}

inline void
StructureTopology::list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist,
        const connection_kind_type kind, std::vector<std::size_t>& out) const
{
    if(dist == 0)
    {
        if(std::find(out.cbegin(), out.cend(), node_idx) != out.cend())
        {
            out.push_back(node_idx);
        }
        return;
    }

    for(auto ik : this->nodes_.at(node_idx).adjacents)
    {
        if(ik.second != kind ||
           std::find(out.cbegin(), out.cend(), ik.first) != out.cend())
        {
            continue;
        }
        out.push_back(ik.first);
        this->list_adjacent_within(ik.first, dist-1, kind, out);
    }
    return;
}

inline void
StructureTopology::construct_chains()
{
    if(this->nodes_.empty()){return;}
    for(auto& node : nodes_)
    {
        node.chain_id = uninitialized();
    }

    std::size_t next_chain_id = 0;
    for(auto& node : nodes_)
    {
//         node.chain_id = uninitialized(); // already set as `uninit`
        for(const auto& edge : node.adjacents)
        {
            if(edge.second == connection_kind_type::contact) {continue;}

            node.chain_id = this->nodes_.at(edge.first).chain_id;
            if(node.chain_id != uninitialized()) {break;}
        }
        if(node.chain_id == uninitialized())
        {
            node.chain_id = next_chain_id++;
        }
    }

    this->num_chains_ = ++next_chain_id;
    return;
}

} // mjolnir
#endif// MJOLNIR_STRUCTURE_TOPOLOGY_H
