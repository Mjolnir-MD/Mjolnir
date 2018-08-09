#ifndef MJOLNIR_MATH_CONSTANTS_HPP
#define MJOLNIR_MATH_CONSTANTS_HPP

namespace mjolnir
{
namespace math
{

template<typename realT>
struct constants
{
    typedef realT real_type;
    static constexpr real_type tolerance = static_cast<real_type>(1e-8);

    static constexpr real_type pi          = static_cast<real_type>(3.14159265358979323846264338);
    static constexpr real_type half_pi     = static_cast<real_type>(1.57079632679489661923132169);
    static constexpr real_type two_pi      = static_cast<real_type>(6.28318530717958647692528677);
    static constexpr real_type one_over_pi = static_cast<real_type>(0.31830988618379067153776753);
};

template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::tolerance;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::pi;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::half_pi;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::two_pi;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::one_over_pi;

} // math
} // mjolnir
#endif// MJOLNIR_MATH_CONSTANTS_HPP
