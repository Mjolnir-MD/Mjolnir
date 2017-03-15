#ifndef MEGINGJORD_SIMD_PACK_PLAIN
#define MEGINGJORD_SIMD_PACK_PLAIN

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT USE megingjord/simd/pack_plain.hpp alone"
#endif

namespace megingjord
{
namespace simd
{

template<>
struct pack<float, 1>
{
    typedef float value_type;
    typedef float type;
    typedef aligned_array<value_type, 1> array_type;
    constexpr static std::size_t size       = 1;
    constexpr static std::size_t align_byte = 4;
};

template<>
struct pack<double, 1>
{
    typedef double value_type;
    typedef double type;
    typedef aligned_array<value_type, 1> array_type;
    constexpr static std::size_t size       = 1;
    constexpr static std::size_t align_byte = 8;
};

struct plain_traits
{
    template<typename T> struct pack_trait{};
    constexpr static std::size_t align_byte = 8;
};

template<>
struct plain_traits::pack_trait<double>
{
    typedef pack<double, 1> type;
};
template<>
struct plain_traits::pack_trait<float>
{
    typedef pack<float, 1> type;
};

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_PACK_AVX */
