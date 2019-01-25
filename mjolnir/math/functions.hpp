#ifndef MJOLNIR_MATH_FUNCTIONS_HPP
#define MJOLNIR_MATH_FUNCTIONS_HPP
#include <cstdint>
#include <cmath>

#if defined(__SSE__) && defined(__GNUC__)
#include <x86intrin.h>
#endif

namespace mjolnir
{
namespace math
{

// std::clamp is defined after C++17
// after C++14, std::min and std::max become constexpr
template<typename T>
inline T clamp(const T x, const T low, const T high) noexcept
{
    return std::min(std::max(x, low), high);
}


// ---------------------------------------------------------------------------
// rsqrt
// ---------------------------------------------------------------------------

template<typename realT>
inline realT rsqrt(realT x) noexcept
{
    return realT(1.0) / std::sqrt(x);
}

#if defined(MJOLNIR_WITH_APPROX)
#  if defined(__SSE__)

template<>
inline float rsqrt<float>(float x) noexcept
{
    return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(x)));
}

#  endif // SSE
#  if defined(__AVX512F__) && defined(__SSE2__)

template<>
inline double rsqrt<double>(double x) noexcept
{
    return _mm_cvtsd_f64(_mm_rsqrt14_sd(_mm_undefined_pd(), _mm_set_sd(x)));
}

#  endif // AVX512F
#endif // MJOLNIR_WITH_APPROX

} // math
} // mjolnir
#endif // MJOLNIR_MATH_FUNCTIONS_HPP
