#ifndef MJOLNIR_GROUP_IGNORATION
#define MJOLNIR_GROUP_IGNORATION

namespace mjolnir
{

template<typename group_idT>
struct GroupIgnorationBase
{
    typedef group_idT group_id_type;
    virtual ~GroupIgnorationBase() = default;
    virtual bool is_ignored(const group_id_type& i, const group_id_type& j
                            ) const noexcept = 0;
};

template<typename group_idT>
struct IgnoreNothing final : public GroupIgnorationBase<group_idT>
{
    typedef GroupIgnorationBase<group_idT> base_type;
    typedef typename base_type::group_id_type group_id_type;
    IgnoreNothing() = default;
    ~IgnoreNothing() override = default;
    bool is_ignored(const group_id_type& i, const group_id_type& j
                    ) const noexcept override
    {
        return false;
    }
};

template<typename group_idT>
struct IgnoreSelf final : public GroupIgnorationBase<group_idT>
{
    typedef GroupIgnorationBase<group_idT> base_type;
    typedef typename base_type::group_id_type group_id_type;
    IgnoreSelf() = default;
    ~IgnoreSelf() override = default;
    bool is_ignored(const group_id_type& i, const group_id_type& j
                    ) const noexcept override
    {
        return i == j;
    }
};

template<typename group_idT>
struct IgnoreOthers final : public GroupIgnorationBase<group_idT>
{
    typedef GroupIgnorationBase<group_idT> base_type;
    typedef typename base_type::group_id_type group_id_type;
    IgnoreOthers() = default;
    ~IgnoreOthers() override = default;
    bool is_ignored(const group_id_type& i, const group_id_type& j
                    ) const noexcept override
    {
        return i != j;
    }
};

} // mjolnir
#endif// MJOLNIR_GROUP_IGNORATION
