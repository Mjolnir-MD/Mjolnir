#ifndef MEGINGJORD_SIMD_OPERATION_ARRAY
#define MEGINGJORD_SIMD_OPERATION_ARRAY

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT USE this header alone"
#endif

#include <immintrin.h>

namespace megingjord
{
namespace simd
{

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S>
operator+(const packed_array<T, N, S>& lhs, const packed_array<T, N, S>& rhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = add_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i], rhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S>
operator-(const packed_array<T, N, S>& lhs, const packed_array<T, N, S>& rhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = sub_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i], rhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S>
operator*(const packed_array<T, N, S>& lhs, const packed_array<T, N, S>& rhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = mul_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i], rhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S>
operator/(const packed_array<T, N, S>& lhs, const packed_array<T, N, S>& rhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = div_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i], rhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S> rcp(const packed_array<T, N, S>& lhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = rcp_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S> rsqrt(const packed_array<T, N, S>& lhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = rsqrt_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S> sqrt(const packed_array<T, N, S>& lhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = sqrt_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S> floor(const packed_array<T, N, S>& lhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = floor_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S> ceil(const packed_array<T, N, S>& lhs)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = ceil_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(lhs[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S>
fmadd(const packed_array<T, N, S>& a, const packed_array<T, N, S>& b,
      const packed_array<T, N, S>& c)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = fmadd_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(a[i], b[i], c[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S>
fmsub(const packed_array<T, N, S>& a, const packed_array<T, N, S>& b,
      const packed_array<T, N, S>& c)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = fmsub_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(a[i], b[i], c[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S>
fnmadd(const packed_array<T, N, S>& a, const packed_array<T, N, S>& b,
       const packed_array<T, N, S>& c)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = fnmadd_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(a[i], b[i], c[i]);
    return retval;
}

template<typename T, std::size_t N, std::size_t S>
inline packed_array<T, N, S>
fnmsub(const packed_array<T, N, S>& a, const packed_array<T, N, S>& b,
       const packed_array<T, N, S>& c)
{
    packed_array<T, N, S> retval;
    for(std::size_t i=0; i<N; ++i)
        retval[i] = fnmsub_impl<typename packed_array<T, N, S>::packed_type
            >::invoke(a[i], b[i], c[i]);
    return retval;
}

// TODO: load(?), store, set(?), broadcast


template<typename T, std::size_t N, std::size_t S>
inline typename packed_array<T, N, S>::aligned_array_type
store(const packed_array<T, N, S>& pa)
{
    if(packed_array<T, N, S>::filled) //TODO do this at compile time
    {
        typename packed_array<T, N, S>::aligned_array_type retval;
        T* ptr = retval.data();
        for(std::size_t i=0; i<packed_array<T, N, S>::packed_size; ++i)
        {
            store_impl<typename packed_array<T, N, S>::packed_type>::invoke(
                    ptr, pa[i]);
            ptr += packed_array<T, N, S>::pack_size;
        }
        return retval;
    }
    else
    {
        typename packed_array<T, N, S>::filled_array_type tmp;
        typename packed_array<T, N, S>::aligned_array_type retval;
        T* ptr = tmp.data();
        for(std::size_t i=0; i<packed_array<T, N, S>::packed_size-1; ++i)
        {
            store_impl<typename packed_array<T, N, S>::packed_type>::invoke(
                    ptr, pa[i]);
            ptr += packed_array<T, N, S>::pack_size;
        }

        for(std::size_t i=0; i<N; ++i)
        {
            retval[i] = tmp[i];
        }
        return retval;
    }
}









} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_OPERATION_ARRAY */
