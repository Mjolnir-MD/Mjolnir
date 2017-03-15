#ifndef MEGINGJORD_SIMD_OPERATION_ARRAY
#define MEGINGJORD_SIMD_OPERATION_ARRAY

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT USE this header alone"
#endif

#include "operation_pack.hpp"
#include "packable_array.hpp"

namespace megingjord
{
namespace simd
{

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
operator+(const packable_array<T, N, S>& lhs, const packable_array<T, N, S>& rhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;

    auto liter = lhs.cpbegin();
    auto riter = rhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            add_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(liter.raw()),
                load_impl<packed_type>::invoke(riter.raw()))
            );
        ++liter; ++riter; ++iter;
    }
    return retval;
}


template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
operator-(const packable_array<T, N, S>& lhs, const packable_array<T, N, S>& rhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;

    auto liter = lhs.cpbegin();
    auto riter = rhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            sub_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(liter.raw()),
                load_impl<packed_type>::invoke(riter.raw()))
            );
        ++liter; ++riter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
operator*(const packable_array<T, N, S>& lhs, const packable_array<T, N, S>& rhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;

    auto liter = lhs.cpbegin();
    auto riter = rhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            mul_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(liter.raw()),
                load_impl<packed_type>::invoke(riter.raw()))
            );
        ++liter; ++riter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
operator*(const packable_array<T, N, S>& lhs, const T rhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto liter = lhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            mul_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(liter.raw()),
                set_impl<packed_type>::invoke(rhs))
            );
        ++liter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
operator*(const T lhs, const packable_array<T, N, S>& rhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto riter = rhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            mul_impl<packed_type>::invoke(
                set_impl<packed_type>::invoke(lhs),
                load_impl<packed_type>::invoke(riter.raw()))
            );
        ++riter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
operator/(const packable_array<T, N, S>& lhs, const packable_array<T, N, S>& rhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;

    auto liter = lhs.cpbegin();
    auto riter = rhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            div_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(liter.raw()),
                load_impl<packed_type>::invoke(riter.raw()))
            );
        ++liter; ++riter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
operator/(const packable_array<T, N, S>& lhs, const T rhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto liter = lhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            div_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(liter.raw()),
                set_impl<packed_type>::invoke(rhs))
            );
        ++liter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
operator/(const T lhs, const packable_array<T, N, S>& rhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto riter = rhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            div_impl<packed_type>::invoke(
                set_impl<packed_type>::invoke(lhs),
                load_impl<packed_type>::invoke(riter.raw()))
            );
        ++riter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S> floor(const packable_array<T, N, S>& lhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto liter = lhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            floor_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(liter.raw()))
            );
        ++liter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S> ceil(const packable_array<T, N, S>& lhs)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto liter = lhs.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            ceil_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(liter.raw()))
            );
        ++liter; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
fmadd(const packable_array<T, N, S>& a, const packable_array<T, N, S>& b,
      const packable_array<T, N, S>& c)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto aiter = a.cpbegin();
    auto biter = b.cpbegin();
    auto citer = c.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            fmadd_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(aiter.raw()),
                load_impl<packed_type>::invoke(biter.raw()),
                load_impl<packed_type>::invoke(citer.raw())
                )
            );
        ++aiter; ++biter; ++citer; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
fmsub(const packable_array<T, N, S>& a, const packable_array<T, N, S>& b,
      const packable_array<T, N, S>& c)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto aiter = a.cpbegin();
    auto biter = b.cpbegin();
    auto citer = c.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            fmsub_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(aiter.raw()),
                load_impl<packed_type>::invoke(biter.raw()),
                load_impl<packed_type>::invoke(citer.raw())
                )
            );
        ++aiter; ++biter; ++citer; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
fnmadd(const packable_array<T, N, S>& a, const packable_array<T, N, S>& b,
       const packable_array<T, N, S>& c)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto aiter = a.cpbegin();
    auto biter = b.cpbegin();
    auto citer = c.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            fnmadd_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(aiter.raw()),
                load_impl<packed_type>::invoke(biter.raw()),
                load_impl<packed_type>::invoke(citer.raw())
                )
            );
        ++aiter; ++biter; ++citer; ++iter;
    }
    return retval;
}

template<typename T, std::size_t N, typename S>
inline packable_array<T, N, S>
fnmsub(const packable_array<T, N, S>& a, const packable_array<T, N, S>& b,
       const packable_array<T, N, S>& c)
{
    typedef typename packable_array<T, N, S>::packed_type packed_type;
    packable_array<T, N, S> retval;
    auto aiter = a.cpbegin();
    auto biter = b.cpbegin();
    auto citer = c.cpbegin();
    auto iter  = retval.pbegin();

    for(std::size_t i=0; i<packable_array<T, N, S>::number_of_packs; ++i)
    {
        store_impl<packed_type>::invoke(iter.raw(),
            fnmsub_impl<packed_type>::invoke(
                load_impl<packed_type>::invoke(aiter.raw()),
                load_impl<packed_type>::invoke(biter.raw()),
                load_impl<packed_type>::invoke(citer.raw())
                )
            );
        ++aiter; ++biter; ++citer; ++iter;
    }
    return retval;
}

// rcp, rsqrt, sqrt ?

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_OPERATION_ARRAY */
