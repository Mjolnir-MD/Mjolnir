#ifndef JARNGREIPR_FORCEFIELD_CONNECTING_INDICES
#define JARNGREIPR_FORCEFIELD_CONNECTING_INDICES
#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

namespace mjolnir
{

struct ConnectingIndices
{
    typedef std::vector<std::size_t> container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

    ConnectionIndices() = default;
    ~ConnectionIndices() = default;

    bool insert(std::size_t i)
    {
        if(indices_.cend() == std::find(indices_.cbegin(), indices_.cend(), i))
        {
            indices_.push_back(i);
            std::sort(indices_.begin(), indices_.end());
        }
        return false;
    }

    bool erase(std::size_t i)
    {
        const auto found = std::find(indices_.cbegin(), indices_.cend(), i);
        if(indices_.cend() != found)
        {
            indices_.erase(found);
            assert(std::is_sorted(indices_.cbegin(), indices_.cend()));
            return true;
        }
        return false;
    }

    void merge(const ConnectingIndices& rhs)
    {
        if(this->empty())
        {
            this->indices_ = rhs;
            return;
        }

        std::vector<std::size_t> tmp(this->size() + rhs.size());
        std::merge(this->cbegin(), this->cend(), rhs.cbegin(), rhs.cend(),
                   tmp.begin());
        this->indices_.clear();
        std::unique_copy(tmp.begin(), tmp.end(),
                         std::back_inserter(this->indices_));
        return;
    }

    std::size_t size() const noexcept {return indices_.size();}
    bool empty() const noexcept {return indices_.size();}
    void clear(){indices_.clear();}

    iterator       begin()        noexcept {return indices_.begin();}
    iterator       end()          noexcept {return indices_.end();}
    const_iterator begin()  const noexcept {return indices_.begin();}
    const_iterator end()    const noexcept {return indices_.end();}
    const_iterator cbegin() const noexcept {return indices_.cbegin();}
    const_iterator cend()   const noexcept {return indices_.cend();}

  private:
    std::vector<std::size_t> indices_;
};


template<typename Map>
Map merge_connection_info(const Map& lhs, const Map& rhs)
{
    static_assert(
            std::is_same<typename Map::value_type, ConnectionIndices>::value,
            "jarngreipr/forcefield/ConnectionInfo.hpp: merge_connection_info: "
            "this function is for map<size_t, connection_indices>.");
    Map result;
    for(const auto& kv : lhs)
    {
        result[kv.first].merge(kv.second);
    }
    for(const auto& kv : rhs)
    {
        result[kv.first].merge(kv.second);
    }
    return result;
}

} // mjolnir
#endif// JARNGREIPR_FORCEFIELD_CONNECTION_INFO
