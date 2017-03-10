#ifndef MEGINGJORD_SIMD_PACK
#define MEGINGJORD_SIMD_PACK
#include <megingjord/util/aligned_array.hpp>
#include <cmath>

namespace megingjord
{
namespace simd
{

template<typename, std::size_t> struct pack;
template<typename, std::size_t, std::size_t> struct packed_array;

template<typename> struct is_simd : public std::false_type{};
template<typename> struct single_type_of;
template<typename T, std::size_t N>
struct single_type_of<pack<T, N>>
{
    typedef T type;
};
template<typename T, std::size_t N, std::size_t A>
struct single_type_of<aligned_array<T, N, A>>
{
    typedef T type;
};
template<typename T, std::size_t N, std::size_t A>
struct single_type_of<packed_array<T, N, A>>
{
    typedef T type;
};

template<typename> struct set_impl;       // from register to register
template<typename> struct load_impl;      // from memory to simd-register(full)
template<typename> struct broadcast_impl; // from memory to simd-register(fewer)
template<typename> struct store_impl;     // from simd-register to memory

template<typename T> struct add_impl
{
    constexpr static inline
    T invoke(T lhs, T rhs) {return lhs + rhs;}
};

template<typename T> struct sub_impl
{
    constexpr static inline
    T invoke(T lhs, T rhs) {return lhs - rhs;}
};

template<typename T> struct mul_impl
{
    constexpr static inline
    T invoke(T lhs, T rhs) {return lhs * rhs;}
};

template<typename T> struct div_impl
{
    constexpr static inline
    T invoke(T lhs, T rhs) {return lhs / rhs;}
};

template<typename T> struct rcp_impl
{
    constexpr static inline
    T invoke(T a) {return 1 / a;}
};

template<typename T> struct rsqrt_impl
{
    static inline
    T invoke(T a) {return 1 / std::sqrt(a);}
};

template<typename T> struct sqrt_impl
{
    static inline
    T invoke(T a) {return std::sqrt(a);}
};

template<typename T> struct floor_impl
{
    static inline
    T invoke(T a) {return std::floor(a);}
};

template<typename T> struct ceil_impl
{
    static inline
    T invoke(T a) {return std::ceil(a);}
};

template<typename T> struct fmadd_impl
{
    constexpr static inline
    T invoke(T a, T b, T c) {return a * b + c;}
};

template<typename T> struct fnmadd_impl
{
    constexpr static inline
    T invoke(T a, T b, T c) {return -a * b + c;}
};

template<typename T> struct fmsub_impl
{
    constexpr static inline
    T invoke(T a, T b, T c) {return a * b - c;}
};

template<typename T> struct fnmsub_impl
{
    constexpr static inline
    T invoke(T a, T b, T c) {return -a * b - c;}
};

} // simd
} // megingjord

#ifdef MJOLNIR_HAVE_AVX
#include <immintrin.h>
#include "avx/pack_avx.hpp"
#include "avx/functor_avx.hpp"
#endif

#include "operation_pack.hpp"
#include "packed_array.hpp"
#include "operation_array.hpp"

#endif /* MEGINGJORD_SIMD_PACK */
