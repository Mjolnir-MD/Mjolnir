#ifndef MJOLNIR_CORE_UNIT_HPP
#define MJOLNIR_CORE_UNIT_HPP
#include <string>

namespace mjolnir
{

namespace unit
{

// -----------------------------------------------------------------------------
// At first the unit system is implemented based on phantom-type. However,
// since mjolnir doesn't have any plan to use such a complex unit system
// inside the engine, we replaced them to simpler one, the constants to convert.
// -----------------------------------------------------------------------------

template<typename realT>
struct constants
{
    using real_type = realT;

    // NIST CODATA in the SI unit
    static constexpr real_type boltzmann_constant()  {return 1.38064852e-23;  } // [J/K]
    static constexpr real_type avogadro_constant()   {return 6.022140857e23;  } // [/mol]
    static constexpr real_type elementary_charge()   {return 1.6021766208e-19;} // [C]
    static constexpr real_type vacuum_permittivity() {return 8.854187817e-12; } // [F/m]
    // [V] = [J/C]; [F/m] = [C/V.m] = [C^2/J.m]

    // ------------------------------------------------------------------------
    // constants to convert units
    // ------------------------------------------------------------------------

    // 1.0 [nm] * nm_to_angstrom = 10.0 [A]
    static constexpr real_type nm_to_angstrom() {return 1e+1;}  // 1 nm = 10 A
    static constexpr real_type angstrom_to_nm() {return 1e-1;}  // 1 A  = 0.1 nm

    static constexpr real_type m_to_angstrom() {return 1e+10;} // 1 m  = 10^+10 a
    static constexpr real_type angstrom_to_m() {return 1e-10;} // 1 A  = 10^-10 A

    static constexpr real_type m_to_nm() {return 1e+9;}  // 1 m  = 10^+9 nm
    static constexpr real_type nm_to_m() {return 1e-9;}  // 1 nm = 1e^-9 m

    static constexpr real_type cal_to_J() {return 4.184;}
    static constexpr real_type J_to_cal() {return 1.0 / cal_to_J();}
};
} // unit

// -----------------------------------------------------------------------------
// List of physical constants. The values may be modified after reading the
// input files by using conversion coefficients. Before that, the values have
// the default values defined in the unit::constants.
// -----------------------------------------------------------------------------

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

    static void set_energy_unit(const std::string& name)
    {
        energy_unit_ = name;
    }
    static void set_length_unit(const std::string& name)
    {
        length_unit_ = name;
    }
    static std::string const& energy_unit() noexcept
    {
        return energy_unit_;
    }
    static std::string const& length_unit() noexcept
    {
        return length_unit_;
    }

    // forget all the units and put the default value. mainly for test codes
    static void reset() noexcept
    {
        kB_   = unit_type::boltzmann_constant();
        NA_   = unit_type::avogadro_constant();
        eps0_ = unit_type::vacuum_permittivity() /
                unit_type::elementary_charge()   /
                unit_type::elementary_charge();

        length_unit_ = "m";
        energy_unit_ = "J";

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

    static std::string energy_unit_;
    static std::string length_unit_;
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
realT constants<realT>::kB_ = unit::constants<realT>::boltzmann_constant();
template<typename realT>
realT constants<realT>::NA_ = unit::constants<realT>::avogadro_constant();
template<typename realT>
realT constants<realT>::eps0_ = (unit::constants<realT>::vacuum_permittivity() /
                                 unit::constants<realT>::elementary_charge())  /
                                 unit::constants<realT>::elementary_charge();

template<typename realT> realT constants<realT>::m_to_length_ = 1.0;
template<typename realT> realT constants<realT>::length_to_m_ = 1.0;
template<typename realT> realT constants<realT>::L_to_volume_ = 1e-3;
template<typename realT> realT constants<realT>::volume_to_L_ = 1e+3;
template<typename realT> std::string constants<realT>::length_unit_ = "J";
template<typename realT> std::string constants<realT>::energy_unit_ = "m";

} // physics

} // mjolnir
#endif// MJOLNIR_CORE_UNIT_HPP
