#ifndef MEGINGJORD_SIMD_FUNCTOR_AVX
#define MEGINGJORD_SIMD_FUNCTOR_AVX
#include <type_traits>

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT USE megingjord/simd/pack_avx.hpp alone"
#endif

#ifndef MJOLNIR_HAVE_AVX
#error "no AVX on the architecture"
#endif

namespace megingjord
{
namespace simd
{

template<> struct is_simd<__m256>  : public std::true_type {};
template<> struct is_simd<__m256d> : public std::true_type {};

template<>
struct add_impl<__m256>
{
    static inline
    __m256 invoke(__m256 x, __m256 y)
    {
        return _mm256_add_ps(x, y);
    }
};

template<>
struct add_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d x, __m256d y)
    {
        return _mm256_add_pd(x, y);
    }
};


template<>
struct sub_impl<__m256>
{
    static inline
    __m256 invoke(__m256 x, __m256 y)
    {
        return _mm256_sub_ps(x, y);
    }
};

template<>
struct sub_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d x, __m256d y)
    {
        return _mm256_sub_pd(x, y);
    }
};

template<>
struct mul_impl<__m256>
{
    static inline
    __m256 invoke(__m256 x, __m256 y)
    {
        return _mm256_mul_ps(x, y);
    }
};

template<>
struct mul_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d x, __m256d y)
    {
        return _mm256_mul_pd(x, y);
    }
};

template<>
struct div_impl<__m256>
{
    static inline
    __m256 invoke(__m256 x, __m256 y)
    {
        return _mm256_div_ps(x, y);
    }
};

template<>
struct div_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d x, __m256d y)
    {
        return _mm256_div_pd(x, y);
    }
};

template<>
struct ceil_impl<__m256>
{
    static inline
    __m256 invoke(__m256 x)
    {
        return _mm256_ceil_ps(x);
    }
};

template<>
struct ceil_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d x)
    {
        return _mm256_ceil_pd(x);
    }
};

template<>
struct floor_impl<__m256>
{
    static inline
    __m256 invoke(__m256 x)
    {
        return _mm256_floor_ps(x);
    }
};

template<>
struct floor_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d x)
    {
        return _mm256_floor_pd(x);
    }
};

// TODO: emulate rcp_pd, rsqrt_pd
template<>
struct rcp_impl<__m256>
{
    static inline
    __m256 invoke(__m256 a)
    {
        return _mm256_rcp_ps(a);
    }
};

template<>
struct rsqrt_impl<__m256>
{
    static inline
    __m256 invoke(__m256 a)
    {
        return _mm256_rsqrt_ps(a);
    }
};

template<>
struct sqrt_impl<__m256>
{
    static inline
    __m256 invoke(__m256 a)
    {
        return _mm256_sqrt_ps(a);
    }
};

template<>
struct sqrt_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d a)
    {
        return _mm256_sqrt_pd(a);
    }
};

template<>
struct set_impl<__m256>
{
    static inline
    __m256 invoke(float e7, float e6, float e5, float e4,
                  float e3, float e2, float e1, float e0)
    {
        return _mm256_set_ps(e7, e6, e5, e4, e3, e2, e1, e0);
    }

    static inline
    __m256 invoke(float x)
    {
        return _mm256_set1_ps(x);
    }
};

template<>
struct set_impl<__m256d>
{
    static inline
    __m256d invoke(double e3, double e2, double e1, double e0)
    {
        return _mm256_set_pd(e3, e2, e1, e0);
    }

    static inline
    __m256d invoke(double x)
    {
        return _mm256_set1_pd(x);
    }
};

// specific version for pack<double, 4>
// XXX arguments are reversed!
inline typename pack<double, 4>::type
set(double x0, double x1, double x2, double x3)
{
    return set_impl<typename pack<double, 4>::type>::invoke(x3, x2, x1, x0);
}

inline typename pack<float, 8>::type
set(float e0, float e1, float e2, float e3,
    float e4, float e5, float e6, float e7)
{
    return set_impl<typename pack<float, 8>::type>::invoke(
            e7, e6, e5, e4, e3, e2, e1, e0);
}

template<>
struct load_impl<__m256>
{
    static inline
    __m256 invoke(const aligned_array<float, 8, 32>& adr)
    {
        return _mm256_load_ps(adr.data());
    }

    static inline
    __m256 invoke(const float* adr)
    {
        return _mm256_load_ps(adr);
    }
};

template<>
struct load_impl<__m256d>
{
    static inline
    __m256d invoke(const aligned_array<double, 4, 32>& adr)
    {
        return _mm256_load_pd(adr.data());
    }

    static inline
    __m256d invoke(const double* adr)
    {
        return _mm256_load_pd(adr);
    }
};


template<>
struct broadcast_impl<__m256>
{
    static inline
    __m256 invoke(const __m128* adr)
    {
        return _mm256_broadcast_ps(adr);
    }

    static inline
    __m256 invoke(const float* adr)
    {
        return _mm256_broadcast_ss(adr);
    }
};

template<>
struct broadcast_impl<__m256d>
{
    static inline
    __m256d invoke(const __m128d* adr)
    {
        return _mm256_broadcast_pd(adr);
    }

    static inline
    __m256d invoke(const double* adr)
    {
        return _mm256_broadcast_sd(adr);
    }
};

template<>
struct store_impl<__m256>
{
    static inline
    void invoke(float* adr, __m256 x)
    {
        _mm256_store_ps(adr, x);
        return;
    }

    static inline
    aligned_array<float, 8, 32> invoke(__m256 x)
    {
        aligned_array<float, 8, 32> retval;
        _mm256_store_ps(retval.data(), x);
        return retval;
    }
};

template<>
struct store_impl<__m256d>
{
    static inline
    void invoke(double* adr, __m256d x)
    {
        _mm256_store_pd(adr, x);
        return;
    }

    static inline
    aligned_array<double, 4, 32> invoke(__m256d x)
    {
        aligned_array<double, 4, 32> retval;
        _mm256_store_pd(retval.data(), x);
        return retval;
    }
};

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_PACK_AVX */
