#ifndef MJOLNIR_PHYSICAL_CONSTANTS_HPP
#define MJOLNIR_PHYSICAL_CONSTANTS_HPP

namespace mjolnir
{
namespace physics
{

// TODO consider units

template<typename realT>
struct constants
{
    typedef realT real_type;
    static real_type kB;
    static real_type NA;
    static real_type e;
    static real_type vacuum_permittivity;
};

template<typename realT>
typename constants<realT>::real_type constants<realT>::kB = 1.986231313e-3;
template<typename realT>
typename constants<realT>::real_type constants<realT>::NA = 6.0221417930e23;
template<typename realT>
typename constants<realT>::real_type constants<realT>::e = 1.60217648740e-19;
template<typename realT>
typename constants<realT>::real_type constants<realT>::vacuum_permittivity = 8.854187817e-12;

} // physics
} // mjolnir
#endif /* MJOLNIR_CONSTANTS */
