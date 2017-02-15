#ifndef MJOLNIR_CONSTANTS
#define MJOLNIR_CONSTANTS

namespace mjolnir
{

template<typename traitsT>
struct constants
{
    typedef typename traitsT::real_type real_type;
    constexpr static real_type tolerance = 1e-8;
    constexpr static real_type pi        = 3.14159265358979323846264338;
};

template<typename traitsT>
class physics
{
    typedef typename traitsT::real_type real_type;

    // constexpr...?
    static real_type kB;
    static real_type NA;
    static real_type e;
    static real_type vacuum_permittivity;
};

template<typename traitsT>
typename physics<traitsT>::real_type physics<traitsT>::kB = 1.986231313e-3;
template<typename traitsT>
typename physics<traitsT>::real_type physics<traitsT>::NA = 6.0221417930e23;
template<typename traitsT>
typename physics<traitsT>::real_type physics<traitsT>::e = 1.60217648740e-19;
template<typename traitsT>
typename physics<traitsT>::real_type physics<traitsT>::vacuum_permittivity = 8.854187817e-12;

}//mjolnir
#endif /* MJOLNIR_CONSTANTS */
