#ifndef JARNGREIPR_FORCEFIELD_CONNECTING_INDICES
#define JARNGREIPR_FORCEFIELD_CONNECTING_INDICES
#include <algorithm>
#include <vector>
#include <map>

namespace mjolnir
{

struct ConnectingIndices
{
    ConnectionIndex() = default;
    ~ConnectionIndex() = default;

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
            return true;
        }
        return false;
    }

    void clear(){indices_.clear();}

  private:
    std::vector<std::size_t> indices_;
};

} // mjolnir
#endif// JARNGREIPR_FORCEFIELD_CONNECTION_INFO
