#ifndef MJOLNIR_MATH_RSQRT
#define MJOLNIR_MATH_RSQRT
#include <cstdint>

/* Here, functions that calculates reverse square-root is defined. *
 * `rsqrt(x);` returns `1. / sqrt(x)`.                             */

#ifdef MJOLNIR_HAVE_SSE
#include <xmmintrin.h>

namespace mjolnir
{

inline float rsqrt(float x) noexcept
{
    float r;
    _mm_store_ss(&r, _mm_rsqrt_ss(_mm_load_ss(&x)));
    return r * (3.0f - x * r * r) * 0.5f;
}

inline double rsqrt(double x) noexcept
{
    const double xhalf = 0.5 * x;
    float f = static_cast<float>(x);
    _mm_store_ss(&f, _mm_rsqrt_ss(_mm_load_ss(&f)));
    double r = static_cast<double>(f);
    r *= (1.5 - xhalf * r * r);
    r *= (1.5 - xhalf * r * r);
    return r * (1.5 - xhalf * r * r);
}

} // mjolnir

#else // cpu do not have SSE operation

namespace mjolnir
{

inline float rsqrt(float x) noexcept
{
    const float xhalf = 0.5f * x;
    std::int32_t i = *(reinterpret_cast<std::int32_t*>(&x));
    i = 0x5f3759df - (i >> 1);
    x = *(reinterpret_cast<float*>(&i));
    x *= (1.5f - xhalf * x * x);
    return x * (1.5f - xhalf * x * x);
}

inline double rsqrt(double x) noexcept
{
    const double xhalf = 0.5 * x;
    std::int64_t i = *(reinterpret_cast<std::int64_t*>(&x));
    i = 0x5fe6ec85e7de30daLL - (i >> 1);
    x = *(reinterpret_cast<double*>(&i));
    x *= (1.5 - xhalf * x * x);
    x *= (1.5 - xhalf * x * x);
    x *= (1.5 - xhalf * x * x);
    x *= (1.5 - xhalf * x * x);
    return x;
}

}// mjolnir

#endif // HAVE_SSE
#endif // MJOLNIR_MATH_RSQRT
