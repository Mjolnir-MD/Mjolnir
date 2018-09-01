#ifndef MJOLNIR_MATH_RSQRT
#define MJOLNIR_MATH_RSQRT
#include <cmath>

#if defined(__SSE__)
#include <x86intrin.h>
#endif

namespace mjolnir
{

template<typename realT>
inline realT rsqrt(realT x) noexcept
{
    return realT(1.0) / std::sqrt(x);
}

#ifdef __SSE__
template<>
inline float rsqrt<float>(float x) noexcept
{
    return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(x)));
}
#endif

#if defined(__AVX512F__) && defined(__SSE2__)
template<>
inline double rsqrt<double>(double x) noexcept
{
    return _mm_cvtsd_f64(_mm_rsqrt14_sd(_mm_set_sd(x), _mm_undefined_pd()));
}
#endif

} // mjolnir
#endif // MJOLNIR_MATH_RSQRT
