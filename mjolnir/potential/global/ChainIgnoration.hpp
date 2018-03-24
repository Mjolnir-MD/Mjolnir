#ifndef MJOLNIR_CHAIN_IGNORATION
#define MJOLNIR_CHAIN_IGNORATION

namespace mjolnir
{

struct IgnoreNothing
{
     IgnoreNothing() = default;
    ~IgnoreNothing() = default;

    template<typename ChainID>
    bool is_ignored(const ChainID& i, const ChainID& j) const noexcept
    {
        return false;
    }
};

struct IgnoreSelf
{
     IgnoreSelf() = default;
    ~IgnoreSelf() = default;

    template<typename ChainID>
    bool is_ignored(const ChainID& i, const ChainID& j) const noexcept
    {
        return i == j;
    }
};

struct IgnoreOthers
{
     IgnoreOthers() = default;
    ~IgnoreOthers() = default;

    template<typename ChainID>
    bool is_ignored(const ChainID& i, const ChainID& j) const noexcept
    {
        return i != j;
    }
};

} // mjolnir
#endif// MJOLNIR_CHAIN_IGNORATION
