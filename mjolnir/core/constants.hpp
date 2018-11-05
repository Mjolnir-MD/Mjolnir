#ifndef MJOLNIR_PHYSICAL_CONSTANTS_HPP
#define MJOLNIR_PHYSICAL_CONSTANTS_HPP
#include <mjolnir/core/Unit.hpp>
#include <string>

namespace mjolnir
{
namespace physics
{

template<typename realT>
struct constants
{
    typedef realT real_type;
    static real_type kB;   // Boltzmann constant
    static real_type NA;   // Avogadro constant
    static real_type eps0; // vacuum permittivity

    // conversion coefficients
    static real_type m_to_length;
    static real_type length_to_m;

    static real_type L_to_volume;
    static real_type volume_to_L;
};

// convert them in read_units() function defined in mjolnir/input/read_units.hpp
// Here, the unit of electric charge is the elementary charge, not coulomb.
// to convert it, the value of vacuum permittivity is devided by squared
// elementary chage.

// XXX: there are 2 reasons why the unit of charge is defined as the elementary
//      charge. First, usualy charge on particle is defined in such a way (+1
//      means there are 2x elementary charge). Second, the range that single
//      precision `float` can represent is too small (1e+/-38) to store the
//      temporary value while the unit conversion. For example, the value of
//      vacuum permittivity become 6e-42 in some unit system and that exceeds
//      the range.

template<typename realT>
realT constants<realT>::kB = unit::constants<realT>::boltzmann_constant;
template<typename realT>
realT constants<realT>::NA = unit::constants<realT>::avogadro_constant;
template<typename realT>
realT constants<realT>::eps0 = (unit::constants<realT>::vacuum_permittivity /
                                unit::constants<realT>::elementary_charge)  /
                                unit::constants<realT>::elementary_charge;

template<typename realT> realT constants<realT>::m_to_length = 1.0;
template<typename realT> realT constants<realT>::length_to_m = 1.0;
template<typename realT> realT constants<realT>::L_to_volume = 1e-3;
template<typename realT> realT constants<realT>::volume_to_L = 1e+3;

} // physics
} // mjolnir
#endif /* MJOLNIR_CONSTANTS */
