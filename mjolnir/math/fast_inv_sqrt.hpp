#ifndef MJOLNIR_FAST_INV_SQRT
#define MJOLNIR_FAST_INV_SQRT
#include <cstdint>

namespace mjolnir
{
 
inline double fast_inv_sqrt(double x)
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

inline float fast_inv_sqrt(float x)
{
    const float xhalf = 0.5f * x;
    std::int32_t i = *(reinterpret_cast<std::int32_t*>(&x));
    i = 0x5f3759df - (i >> 1);
    x = *(reinterpret_cast<float*>(&i));
    x *= (1.5f - xhalf * x * x);
    return x * (1.5f - xhalf * x * x);
}
 
 
}// mjolnir

#endif /* MJOLNIR_FAST_INV_SQRT */
