#ifndef MEGINGJORD_SIMD_OPERATION
#define MEGINGJORD_SIMD_OPERATION

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT USE this header alone"
#endif

#include <immintrin.h>

namespace megingjord
{
namespace simd
{

// NOTE: consider universal reference + perfect forwarding v.s. copy
//       probably compiler can do the kind of optimization...?
//
// EXAMPLE:
// template<typename T, typename std::enable_if<
//     is_simd<typename std::decay<T>::type>::value, std::nullptr_t>::type = nullptr>
// inline T operator+(T&& lhs, T&& rhs)
// {
//     return add_impl<T>::invoke(std::forward<T>(lhs), std::forward<T>(rhs));
// }

template<typename T, typename std::enable_if<
    is_simd<T>::value, std::nullptr_t>::type = nullptr>
inline T operator+(T lhs, T rhs)
{
    return add_impl<T>::invoke(lhs, rhs);
}

template<typename T, typename std::enable_if<
    is_simd<T>::value, std::nullptr_t>::type = nullptr>
inline T operator-(T lhs, T rhs)
{
    return sub_impl<T>::invoke(lhs, rhs);
}

template<typename T, typename std::enable_if<
    is_simd<T>::value, std::nullptr_t>::type = nullptr>
inline T operator*(T lhs, T rhs)
{
    return mul_impl<T>::invoke(lhs, rhs);
}

template<typename T, typename std::enable_if<
    is_simd<T>::value, std::nullptr_t>::type = nullptr>
inline T operator/(T lhs, T rhs)
{
    return div_impl<T>::invoke(lhs, rhs);
}

template<typename T>
inline T rcp(T x)
{
    return rcp_impl<T>::invoke(x);
}

template<typename T>
inline T rsqrt(T x)
{
    return rsqrt_impl<T>::invoke(x);
}

template<typename T>
inline T sqrt(T x)
{
    return sqrt_impl<T>::invoke(x);
}

template<typename T>
inline T floor(T x)
{
    return floor_impl<T>::invoke(x);
}

template<typename T>
inline T ceil(T x)
{
    return ceil_impl<T>::invoke(x);
}

template<typename T>
inline T fmadd(T x)
{
    return fmadd_impl<T>::invoke(x);
}

template<typename T>
inline T fnmadd(T x)
{
    return fnmadd_impl<T>::invoke(x);
}

template<typename T>
inline T fmsub(T x)
{
    return fmsub_impl<T>::invoke(x);
}

template<typename T>
inline T fnmsub(T x)
{
    return fnmsub_impl<T>::invoke(x);
}

template<typename T, std::size_t N, std::size_t A>
inline typename pack<T, N>::type load(const aligned_array<T, N, A>& xs)
{
    return load_impl<typename pack<T, N>::type>::invoke(xs);
}

template<typename T, std::size_t N>
inline typename pack<T, N>::type broadcast(const T* x)
{
    return broadcast_impl<typename pack<T, N>::type>::invoke(x);
}

template<typename T, std::size_t N, std::size_t A>
inline aligned_array<T, N, A> store(const typename pack<T, N>::type x)
{
    return store_impl<typename pack<T, N>::type>::invoke(x);
}

template<typename T, std::size_t N, std::size_t A>
inline void store(aligned_array<T, N, A>& dst, const typename pack<T, N>::type x)
{
    return store_impl<typename pack<T, N>::type>::invoke(dst.data(), x);
}

template<typename T, std::size_t N>
inline typename pack<T, N>::type set(T x)
{
    return set_impl<typename pack<T, N>::type>::invoke(x);
}

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_OPERATION */
