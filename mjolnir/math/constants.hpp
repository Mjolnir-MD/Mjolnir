#ifndef MJOLNIR_MATH_CONSTANTS_HPP
#define MJOLNIR_MATH_CONSTANTS_HPP

namespace mjolnir
{
namespace math
{

// There are 2 objective to have this class. One is to manage special values
// (e.g. pi) and widely-used variants of them (e.g. 2pi, 1/pi) to reduce
// calculation cost (especially, M_PI is not guaranteed to exist by standard).
// Another objective is to deal with templatized codes. Literal `1.0` becomes
// `double` (you need to write `1.0f` for float literal values). It causes
// template specialization/overload resolution ambiguity when instanciated with
// `float` values. to avoid it, use `constants<real>::two` instead of `2.0`.
template<typename realT>
struct constants
{
    typedef realT real_type;
    static constexpr real_type tolerance = static_cast<real_type>(1e-8);

    static constexpr real_type zero           = static_cast<real_type>(0.0);
    static constexpr real_type one            = static_cast<real_type>(1.0);
    static constexpr real_type two            = static_cast<real_type>(2.0);
    static constexpr real_type three          = static_cast<real_type>(3.0);
    static constexpr real_type four           = static_cast<real_type>(4.0);
    static constexpr real_type five           = static_cast<real_type>(5.0);
    static constexpr real_type six            = static_cast<real_type>(6.0);
    static constexpr real_type seven          = static_cast<real_type>(7.0);
    static constexpr real_type eight          = static_cast<real_type>(8.0);
    static constexpr real_type nine           = static_cast<real_type>(9.0);
    static constexpr real_type ten            = static_cast<real_type>(10.0);

    static constexpr real_type minus_one      = static_cast<real_type>(-1.0);
    static constexpr real_type half           = static_cast<real_type>(0.5);
    static constexpr real_type one_third      = static_cast<real_type>(0.33333333333333333333333333);
    static constexpr real_type two_thirds     = static_cast<real_type>(0.66666666666666666666666667);
    static constexpr real_type quarter        = static_cast<real_type>(0.25);
    static constexpr real_type three_quarters = static_cast<real_type>(0.75);

    static constexpr real_type pi          = static_cast<real_type>(3.14159265358979323846264338);
    static constexpr real_type half_pi     = static_cast<real_type>(1.57079632679489661923132169);
    static constexpr real_type two_pi      = static_cast<real_type>(6.28318530717958647692528677);
    static constexpr real_type one_over_pi = static_cast<real_type>(0.31830988618379067153776753);
};

template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::tolerance;

template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::one;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::two;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::three;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::four;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::five;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::six;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::seven;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::eight;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::nine;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::ten;

template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::minus_one;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::half;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::one_third;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::two_thirds;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::quarter;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::three_quarters;

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
