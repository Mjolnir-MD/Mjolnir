#ifndef MJOLNIR_POTENTIAL_GLOBAL_IGNORE_GROUP_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_IGNORE_GROUP_HPP
#include <mjolnir/util/throw_exception.hpp>
#include <string>
#include <memory>

// These classes are to ignore inter- or intra-group interaction.
//
// These would be contained in a global potential class and ExclusionList
// uses them indirectly.

namespace mjolnir
{

template<typename GroupID>
struct IgnoreGroupBase
{
    using group_id_type = GroupID;

    IgnoreGroupBase()          = default;
    virtual ~IgnoreGroupBase() = default;

    virtual bool is_ignored(
        const group_id_type& i, const group_id_type& j) const noexcept = 0;

    virtual const char* name() const noexcept = 0;
};

// inter group interaction
template<typename GroupID>
struct IgnoreInterGroup final : public IgnoreGroupBase<GroupID>
{
    using base_type     = IgnoreGroupBase<GroupID>;
    using group_id_type = typename base_type::group_id_type;

    IgnoreInterGroup() = default;
    ~IgnoreInterGroup() noexcept override = default;

    bool is_ignored(
        const group_id_type& i, const group_id_type& j) const noexcept override
    {
        // is inter-group interaction
        return i != j;
    }
    const char* name() const noexcept {return "Inter";}
};

// intra group interaction
template<typename GroupID>
struct IgnoreIntraGroup final : public IgnoreGroupBase<GroupID>
{
    using base_type     = IgnoreGroupBase<GroupID>;
    using group_id_type = typename base_type::group_id_type;

    IgnoreIntraGroup() = default;
    ~IgnoreIntraGroup() noexcept override = default;

    bool is_ignored(
        const group_id_type& i, const group_id_type& j) const noexcept override
    {
        // is intra-group interaction
        return i == j;
    }
    const char* name() const noexcept {return "Intra";}
};

// count all the interactions
template<typename GroupID>
struct IgnoreNoGroup final : public IgnoreGroupBase<GroupID>
{
    using base_type     = IgnoreGroupBase<GroupID>;
    using group_id_type = typename base_type::group_id_type;

    IgnoreNoGroup() = default;
    ~IgnoreNoGroup() noexcept override = default;

    bool is_ignored(
        const group_id_type&, const group_id_type&) const noexcept override
    {
        return false;
    }
    const char* name() const noexcept {return "Nothing";}
};

template<typename GroupID>
struct IgnoreGroup
{
    using group_id_type = GroupID;

    explicit IgnoreGroup(const std::string& name): ignore_group_(nullptr)
    {
        this->reset(name);
    }

    template<typename IgnoreSomething>
    explicit IgnoreGroup(std::unique_ptr<IgnoreSomething>&& ptr)
        : ignore_group_(std::move(ptr))
    {}
    ~IgnoreGroup() = default;

    IgnoreGroup(IgnoreGroup const& rhs) {this->reset(rhs.name());}
    IgnoreGroup(IgnoreGroup&&      rhs) = default;
    IgnoreGroup& operator=(IgnoreGroup const& rhs) {this->reset(rhs.name()); return *this;}
    IgnoreGroup& operator=(IgnoreGroup&&      rhs) = default;

    bool is_ignored(const GroupID& i, const GroupID& j) const noexcept
    {
        return ignore_group_->is_ignored(i, j);
    }
    const char* name() const noexcept {return ignore_group_->name();}

    void reset(const std::string& name)
    {
        if(name == "Nothing")
        {
            ignore_group_.reset(new IgnoreNoGroup<GroupID>);
        }
        else if(name == "Self" || name == "Intra")
        {
            ignore_group_.reset(new IgnoreIntraGroup<GroupID>);
        }
        else if(name == "Others" || name == "Inter")
        {
            ignore_group_.reset(new IgnoreInterGroup<GroupID>);
        }
        else
        {
            throw_exception<std::invalid_argument>("IgnoreGroup: ",
                "unknown name appreaed: `" + name +
                "`. Only `Nothing`, `Self|Intra`, or `Others|Intra` are allowed.");
        }
    }

  private:

    std::unique_ptr<IgnoreGroupBase<GroupID>> ignore_group_;
};

} // mjolnir
#endif// MJOLNIR_POTENTIAL_GLOBAL_IGNORE_GROUP_HPP
