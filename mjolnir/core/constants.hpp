#ifndef MJOLNIR_PHYSICAL_CONSTANTS_HPP
#define MJOLNIR_PHYSICAL_CONSTANTS_HPP
#include <mjolnir/core/Unit.hpp>

namespace mjolnir
{
namespace physics
{

// TODO consider units

template<typename realT>
struct constants
{
    typedef realT real_type;
    static real_type kB;   // Boltzmann constant
    static real_type NA;   // Avogadro constant
    static real_type e;    // elementary charge
    static real_type eps0; // vacuum permittivity
};

// initialize with SI units. convert it in read_units() function defined in
// mjolnir/input/read_units.hpp

template<typename realT>
realT constants<realT>::kB = unit::constants<realT>::boltzmann_constant;
template<typename realT>
realT constants<realT>::NA = unit::constants<realT>::avogadro_constant;
template<typename realT>
realT constants<realT>::e  = unit::constants<realT>::elementary_charge;
template<typename realT>
realT constants<realT>::eps0 = unit::constants<realT>::vacuum_permittivity;

} // physics
} // mjolnir
#endif /* MJOLNIR_CONSTANTS */
