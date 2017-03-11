#ifndef MEGINGJORD_SIMD_PACK_AVX
#define MEGINGJORD_SIMD_PACK_AVX

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT USE megingjord/simd/pack_avx.hpp alone"
#endif

#ifndef MJOLNIR_HAVE_AVX
#error "mjolnir detects no AVX on the architecture"
#endif

namespace megingjord
{
namespace simd
{

template<>
struct pack<float, 8>
{
    typedef float  value_type;
    typedef __m256 type;
    typedef aligned_array<value_type, 8> array_type;
    constexpr static std::size_t size       = 8;
    constexpr static std::size_t align_byte = 32;
};

template<>
struct pack<double, 4>
{
    typedef double  value_type;
    typedef __m256d type;
    typedef aligned_array<value_type, 4> array_type;
    constexpr static std::size_t size       = 4;
    constexpr static std::size_t align_byte = 32;
};

struct avx_traits
{
    template<typename T> struct pack_trait{};
    constexpr static std::size_t align_byte = 32;
};

template<>
struct avx_traits::pack_trait<double>
{
    typedef pack<double, 4> type;
};
template<>
struct avx_traits::pack_trait<float>
{
    typedef pack<float, 8> type;
};

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_PACK_AVX */
