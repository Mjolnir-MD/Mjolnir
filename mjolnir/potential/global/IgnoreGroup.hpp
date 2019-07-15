#ifndef MJOLNIR_POTENTIAL_GLOBAL_IGNORE_GROUP_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_IGNORE_GROUP_HPP
#include <mjolnir/core/Topology.hpp>
#include <utility>
#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

// These classes are to ignore inter- or intra-group interaction.
//
// These would be contained in a global potential class and ExclusionList
// uses them indirectly.

namespace mjolnir
{

template<typename GroupID>
class IgnoreGroup
{
  public:

    using group_id_type = GroupID;

    IgnoreGroup(std::map<GroupID, std::vector<GroupID>> ignore)
        : ignores_(std::move(ignore))
    {
        // to make it bidirectional
        for(const auto& kv : this->ignores_)
        {
            for(const auto& ignored : kv.second)
            {
                if(std::find(this->ignores_.at(ignored).begin(),
                             this->ignores_.at(ignored).end(), kv.first) ==
                             this->ignores_.at(ignored).end())
                {
                    this->ignores_.at(ignored).push_back(kv.first);
                }
            }
        }

        for(const auto& kv : this->ignores_)
        {
            for(const auto& ignored : kv.second)
            {
                // should be found
                assert(std::find(this->ignores_.at(ignored).begin(),
                                 this->ignores_.at(ignored).end(), kv.first) !=
                                 this->ignores_.at(ignored).end());
            }
        }
    }

    ~IgnoreGroup() = default;
    IgnoreGroup(IgnoreGroup const& rhs) = default;
    IgnoreGroup(IgnoreGroup&&      rhs) = default;
    IgnoreGroup& operator=(IgnoreGroup const& rhs) = default;
    IgnoreGroup& operator=(IgnoreGroup&&      rhs) = default;

    bool is_ignored(const GroupID& i, const GroupID& j) const noexcept
    {
        const auto ignore_list_of_i = ignores_.find(i);
        if(ignore_list_of_i == ignores_.end())
        {
            return false;
        }

        const auto& ignore_list = ignore_list_of_i->second;

        // ignore if j is found in i
        return std::find(ignore_list.begin(), ignore_list.end(), j) !=
               ignore_list.end();
    }

    // only for testing purpose
    std::map<GroupID, std::vector<GroupID>> const& ignores() const noexcept
    {
        return this->ignores_;
    }

  private:

    std::map<GroupID, std::vector<GroupID>> ignores_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class IgnoreGroup<Topology::group_id_type>;
#endif

} // mjolnir
#endif// MJOLNIR_POTENTIAL_GLOBAL_IGNORE_GROUP_HPP
