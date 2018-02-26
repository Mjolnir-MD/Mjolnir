#ifndef MJOLNIR_CONSTANTS
#define MJOLNIR_CONSTANTS

namespace mjolnir
{

template<typename realT>
struct constants
{
    typedef realT real_type;
    constexpr static real_type tolerance = 1e-8;
    constexpr static real_type pi        = 3.14159265358979323846264338;
};

template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::tolerance;
template<typename realT>
constexpr typename constants<realT>::real_type constants<realT>::pi;


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
