#ifndef MJOLNIR_CONSTANTS
#define MJOLNIR_CONSTANTS

namespace mjolnir
{

// there are 2 objective to have this class. one is to manage special values
// (e.g. pi) and widely-used variants of them (e.g. 2pi, 1/pi).
// another objective is to deal with templatized codes. literal 1.0 becomes
// double. it causes ambiguous template specialization/overload resolution.
// to avoid it, use `constants<real>::two` instead of `2.0`.
template<typename realT>
struct constants
{
    typedef realT real_type;
    static constexpr real_type tolerance = static_cast<real_type>(1e-8);

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


template<typename realT>
struct physics
{
    typedef realT real_type;
    static real_type kB;
    static real_type NA;
    static real_type e;
    static real_type vacuum_permittivity;
};

template<typename realT>
typename physics<realT>::real_type physics<realT>::kB = 1.986231313e-3;
template<typename realT>
typename physics<realT>::real_type physics<realT>::NA = 6.0221417930e23;
template<typename realT>
typename physics<realT>::real_type physics<realT>::e = 1.60217648740e-19;
template<typename realT>
typename physics<realT>::real_type physics<realT>::vacuum_permittivity = 8.854187817e-12;

}//mjolnir
#endif /* MJOLNIR_CONSTANTS */
