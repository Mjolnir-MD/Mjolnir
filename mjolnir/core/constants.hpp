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
    using real_type = realT;
    using unit_type = unit::constants<real_type>;

    // access as constant
    static real_type kB()   noexcept {return kB_;}   // Boltzmann constant
    static real_type NA()   noexcept {return NA_;}   // Avogadro constant
    static real_type eps0() noexcept {return eps0_;} // vacuum permittivity

    static real_type m_to_length() noexcept {return m_to_length_;}
    static real_type length_to_m() noexcept {return length_to_m_;}

    static real_type L_to_volume() noexcept {return L_to_volume_;}
    static real_type volume_to_L() noexcept {return volume_to_L_;}

    // setter to protect constants from modification
    static void set_kB  (const real_type v) noexcept {kB_   = v;}
    static void set_NA  (const real_type v) noexcept {NA_   = v;}
    static void set_eps0(const real_type v) noexcept {eps0_ = v;}
    static void set_m_to_length(const real_type v) noexcept {m_to_length_ = v;}
    static void set_length_to_m(const real_type v) noexcept {length_to_m_ = v;}
    static void set_L_to_volume(const real_type v) noexcept {L_to_volume_ = v;}
    static void set_volume_to_L(const real_type v) noexcept {volume_to_L_ = v;}

    // forget all the units and put the default value. mainly for test codes
    static void reset() noexcept
    {
        kB_   = unit_type::boltzmann_constant;
        NA_   = unit_type::avogadro_constant;
        eps0_ = unit_type::vacuum_permittivity /
                unit_type::elementary_charge   /
                unit_type::elementary_charge;

        m_to_length_ = 1.0;
        length_to_m_ = 1.0;
        L_to_volume_ = 1e-3;
        volume_to_L_ = 1e+3;
    }

    // the code in read_units() can be put here. which is better?

  private:

    static real_type kB_;
    static real_type NA_;
    static real_type eps0_;

    // conversion coefficients
    static real_type m_to_length_;
    static real_type length_to_m_;

    static real_type L_to_volume_;
    static real_type volume_to_L_;
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
realT constants<realT>::kB_ = unit::constants<realT>::boltzmann_constant;
template<typename realT>
realT constants<realT>::NA_ = unit::constants<realT>::avogadro_constant;
template<typename realT>
realT constants<realT>::eps0_ = (unit::constants<realT>::vacuum_permittivity /
                                 unit::constants<realT>::elementary_charge)  /
                                 unit::constants<realT>::elementary_charge;

template<typename realT> realT constants<realT>::m_to_length_ = 1.0;
template<typename realT> realT constants<realT>::length_to_m_ = 1.0;
template<typename realT> realT constants<realT>::L_to_volume_ = 1e-3;
template<typename realT> realT constants<realT>::volume_to_L_ = 1e+3;

} // physics
} // mjolnir
#endif /* MJOLNIR_CONSTANTS */
