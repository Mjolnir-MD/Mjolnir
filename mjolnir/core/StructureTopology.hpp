#ifndef MJOLNIR_STRUCTURE_TOPOLOGY_H
#define MJOLNIR_STRUCTURE_TOPOLOGY_H
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <string>

namespace mjolnir
{

class StructureTopology
{
  public:

    enum connection_kind
    {
        bond,   // bond interaction
        contact // other kinds of contacts
    };

    // using adjacency-list
    struct node
    {
        std::vector<std::pair<std::size_t, connection_kind>> adjacents;
    };

  public:
    StructureTopology()  = default;
    ~StructureTopology() = default;
    StructureTopology(const StructureTopology&) = default;
    StructureTopology(StructureTopology&&)      = default;
    StructureTopology& operator=(const StructureTopology&) = default;
    StructureTopology& operator=(StructureTopology&&)      = default;

    StructureTopology(const std::size_t N): nodes_(N){};

    void resize(const std::size_t N) {nodes_.resize(N); return;}

    void add_connection(const std::size_t i, const std::size_t j,
                        const connection_kind kind);

    void erase_connection(const std::size_t i, const std::size_t j);
    void erase_connection(const std::size_t i, const std::size_t j,
                          const connection_kind kind);

    bool has_connection(const std::size_t i, const std::size_t j);
    bool has_connection(const std::size_t i, const std::size_t j,
                        const connection_kind kind);

    std::vector<std::size_t>
    list_adjacent_within(const std::size_t node_idx, const std::size_t dist) const;
    std::vector<std::size_t>
    list_adjacent_within(const std::size_t node_idx, const std::size_t dist,
                         const connection_kind kind) const;

  private:
    std::vector<node> nodes_;
};

inline void StructureTopology::add_connection(const std::size_t i, const std::size_t j,
        const StructureTopology::connection_kind kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw std::out_of_range(
            std::string("mjolnir::StructureTopology::add_connection: size of nodes = ") +
            std::to_string(nodes_.size()) + std::string(", i = ") +
            std::to_string(i) + std::string(", j = ") + std::to_string(j));
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

inline void StructureTopology::erase_connection(const std::size_t i, const std::size_t j)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw std::out_of_range(
            std::string("mjolnir::StructureTopology::erase_connection: size of nodes = ") +
            std::to_string(nodes_.size()) + std::string(", i = ") +
            std::to_string(i) + std::string(", j = ") + std::to_string(j));
    }

    {
        auto& adjs = nodes_[i].adjacents;
        const auto found = std::find_if(adjs.cbegin(), adjs.cend(),
            [j](const std::pair<std::size_t, connection_kind>& x){
                return x.first == j;
            });
        if(found != adjs.cend())
        {
            adjs.erase(found);
        }
    }
    {
        auto& adjs = nodes_[j].adjacents;
        const auto found = std::find_if(adjs.cbegin(), adjs.cend(),
            [i](const std::pair<std::size_t, connection_kind>& x){
                return x.first == i;
            });
        if(found != adjs.cend())
        {
            adjs.erase(found);
        }
    }
    return ;
}


inline void StructureTopology::erase_connection(const std::size_t i, const std::size_t j,
        const StructureTopology::connection_kind kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw std::out_of_range(
            std::string("mjolnir::StructureTopology::erase_connection: size of nodes = ") +
            std::to_string(nodes_.size()) + std::string(", i = ") +
            std::to_string(i) + std::string(", j = ") + std::to_string(j));
    }

    {
        const auto elem = std::make_pair(j, kind);
        auto& adjs = nodes_[i].adjacents;
        const auto found = std::find(adjs.cbegin(), adjs.cend(), elem);
        if(found != adjs.cend())
        {
            adjs.erase(found);
        }
    }
    {
        const auto elem = std::make_pair(i, kind);
        auto& adjs = nodes_[j].adjacents;
        const auto found = std::find(adjs.cbegin(), adjs.cend(), elem);
        if(found != adjs.cend())
        {
            adjs.erase(found);
        }
    }
    return ;
}

inline bool StructureTopology::has_connection(const std::size_t i, const std::size_t j)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw std::out_of_range(
            std::string("mjolnir::StructureTopology::has_connection: size of nodes = ") +
            std::to_string(nodes_.size()) + std::string(", i = ") +
            std::to_string(i) + std::string(", j = ") + std::to_string(j));
    }

    auto& adjs = nodes_[i].adjacents;
    const auto found = std::find_if(adjs.cbegin(), adjs.cend(),
        [j](const std::pair<std::size_t, connection_kind>& x){
            return x.first == j;
        });

    return found != adjs.cend();
}

inline bool StructureTopology::has_connection(const std::size_t i, const std::size_t j,
                    const connection_kind kind)
{
    if(nodes_.size() <= std::max(i, j))
    {
        throw std::out_of_range(
            std::string("mjolnir::StructureTopology::has_connection: size of nodes = ") +
            std::to_string(nodes_.size()) + std::string(", i = ") +
            std::to_string(i) + std::string(", j = ") + std::to_string(j));
    }

    auto& adjs = nodes_[i].adjacents;
    const auto found = std::find(adjs.cbegin(), adjs.cend(),
                                 std::make_pair(j, kind));
    return found != adjs.cend();
}

inline std::vector<std::size_t>
StructureTopology::list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist) const
{
    if(dist == 0) {return std::vector<std::size_t>{};}

    std::vector<std::size_t> retval;
    for(auto ik : this->nodes_.at(node_idx).adjacents)
    {
        retval.push_back(ik.first);
        if(dist > 1)
        {
            const auto tmp = this->list_adjacent_within(ik.first, dist - 1);
            std::copy(tmp.begin(), tmp.end(), std::back_inserter(retval));
        }
    }
    std::sort(retval.begin(), retval.end());
    const auto new_end = std::unique(retval.begin(), retval.end());
    retval.erase(new_end, retval.end());
    return retval;
}

inline std::vector<std::size_t>
StructureTopology::list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist,
        const connection_kind kind) const
{
    if(dist == 0) {return std::vector<std::size_t>{};}

    std::vector<std::size_t> retval;
    for(auto ik : this->nodes_.at(node_idx).adjacents)
    {
        if(ik.second == kind)
        {
            continue;
        }
        retval.push_back(ik.first);
        if(dist > 1)
        {
            const auto tmp = this->list_adjacent_within(ik.first, dist - 1, kind);
            std::copy(tmp.begin(), tmp.end(), std::back_inserter(retval));
        }
    }
    std::sort(retval.begin(), retval.end());
    const auto new_end = std::unique(retval.begin(), retval.end());
    retval.erase(new_end, retval.end());
    return retval;
}

} // mjolnir
#endif// MJOLNIR_STRUCTURE_TOPOLOGY_H
