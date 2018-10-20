#ifndef MJOLNIR_CORE_UNIT_HPP
#define MJOLNIR_CORE_UNIT_HPP

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

    // 1.0 [nm] * nm_to_angstrom = 10.0 [A]
    static constexpr real_type nm_to_angstrom = 1e+1;  // 1 nm = 10 A
    static constexpr real_type angstrom_to_nm = 1e-1;  // 1 A  = 0.1 nm

    static constexpr real_type m_to_angstrom  = 1e+10; // 1 m  = 10^+10 A
    static constexpr real_type angstrom_to_m  = 1e-10; // 1 A  = 10^-10 A

    static constexpr real_type m_to_nm        = 1e+9;  // 1 m  = 10^+9 nm
    static constexpr real_type nm_to_m        = 1e-9;  // 1 nm = 1e^-9 m

    static constexpr real_type cal_to_J       = 4.1868; // the IT calorie
    static constexpr real_type J_to_cal       = 1.0 / kcal_to_kJ;

    static constexpr real_type J_to_eV        = elementary_charge; // J = CxV
    static constexpr real_type eV_to_J        = 1.0 / J_to_eV;

    static constexpr real_type cal_to_eV      = cal_to_J * J_to_eV;
    static constexpr real_type eV_to_cal      =  eV_to_J * J_to_cal;

    // NIST CODATA in the SI unit
    static constexpr real_type boltzmann_constant  = 1.38064852e-23;   // [J/K]
    static constexpr real_type avogadro_constant   = 6.022140857e23;   // [/mol]
    static constexpr real_type elementary_charge   = 1.6021766208e-19; // [C]
    static constexpr real_type vacuum_permittivity = 8.854187817e-12;  // [F/m]
    // [V] = [J/C]; [F/m] = [C/V.m] = [C^2/J.m]
};

template<typename realT> constexpr realT constants<realT>::m_to_nm;
template<typename realT> constexpr realT constants<realT>::nm_to_m;

template<typename realT> constexpr realT constants<realT>::m_to_angstrom;
template<typename realT> constexpr realT constants<realT>::angstrom_to_m;

template<typename realT> constexpr realT constants<realT>::nm_to_angstrom;
template<typename realT> constexpr realT constants<realT>::angstrom_to_nm;

template<typename realT> constexpr realT constants<realT>::cal_to_J;
template<typename realT> constexpr realT constants<realT>::J_to_cal;

template<typename realT> constexpr realT constants<realT>::J_to_eV;
template<typename realT> constexpr realT constants<realT>::eV_to_J;

template<typename realT> constexpr realT constants<realT>::cal_to_eV;
template<typename realT> constexpr realT constants<realT>::eV_to_cal;

template<typename realT> constexpr realT constants<realT>::boltzmann_constant;
template<typename realT> constexpr realT constants<realT>::avogadro_constant;
template<typename realT> constexpr realT constants<realT>::elementary_charge;
template<typename realT> constexpr realT constants<realT>::vacuum_permittivity;

} // unit
} // mjolnir
#endif// MJOLNIR_CORE_UNIT_HPP
