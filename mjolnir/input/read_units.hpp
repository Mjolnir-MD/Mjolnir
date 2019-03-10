#ifndef MJOLNIR_READ_UNIT_SYSTEM_HPP
#define MJOLNIR_READ_UNIT_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_simulator.hpp>
#include <cmath>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<SimulatorBase> read_units(const toml::table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;
    using phys_type = physics::constants<real_type>;
    using math_type = math::constants<real_type>;
    using unit_type = unit::constants<real_type>;

    const auto& units  = toml::find<toml::value>(data,  "units");
    const auto& energy = toml::find<std::string>(units, "energy");
    const auto& length = toml::find<std::string>(units, "length");
    MJOLNIR_LOG_NOTICE("energy unit is [", energy, ']');
    MJOLNIR_LOG_NOTICE("length unit is [", length, ']');

    if(energy == "kcal/mol")
    {
        // kB [J/K] -> [kcal/mol/K] by * (J to cal) * 1e-3 * (/mol)
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal / 1000.0) *
                          unit_type::avogadro_constant);

        // eps0 [F/m] == [C^2/J/m] -> [C^2/(kcal/mol)/m]
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal) /
                            unit_type::avogadro_constant);
    }
    else if(energy == "kJ/mol")
    {
        // kB [J/K] -> [kJ/mol/K]
        phys_type::set_kB(phys_type::kB() * 1e-3 * unit_type::avogadro_constant);
        // eps0 [F/m] == [C^2/J/m] -> [C^2/kJ/mol/m]
        phys_type::set_eps0(phys_type::eps0() * 1e+3 /
                            unit_type::avogadro_constant);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_units: unknown unit for energy: `",
            toml::find(units, "energy"), "here", {
            "expected value is one of the following.",
            "- \"kcal/mol\"",
            "- \"kJ/mol\""
            }));
    }

    // until here, SI `m` are used as length unit.

    // 1 [mol/m^3] = 1e-3  [mol/L] = 1e-27 [mol/nm^3] = 1e-30 [mol/A^3]
    //               1     [mol/L] = 1e-24 [mol/nm^3] = 1e-27 [mol/A^3]
    if(length == "angstrom" || length == u8"Å")
    {
        // eps0 [C^2/Energy/m] -> [C^2/Energy/Angstrom]
        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom);

        // 1 m = 10^9 nm, 1 nm = 10^-9 m
        phys_type::set_m_to_length(unit_type::m_to_angstrom);
        phys_type::set_length_to_m(unit_type::angstrom_to_m);

        // 1 [L] = 1e-3 [m^3] = 1e+24 [nm^3]; 1 [m^3] = 1e27 [nm^3]
        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom, 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m, 3));

        MJOLNIR_LOG_INFO("1 [m] = ", phys_type::m_to_length(), " [", length, "]");
        MJOLNIR_LOG_INFO("1 [", length, "] = ", phys_type::length_to_m(), " [m]");

        MJOLNIR_LOG_INFO("1 [L] = ", phys_type::L_to_volume(), " [", length, "^3]");
        MJOLNIR_LOG_INFO("1 [", length, "^3] = ", phys_type::volume_to_L(), " [L]");
    }
    else if(length == "nm")
    {
        // eps0 [C^2/Energy/m] -> [C^2/Energy/nm]
        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_nm);

        // 1 m = 10^9 nm, 1 nm = 10^-9 m
        phys_type::set_m_to_length(unit_type::m_to_nm);
        phys_type::set_length_to_m(unit_type::nm_to_m);

        // 1 [L] = 1e-3 [m^3] = 1e+24 [nm^3]; 1 [m^3] = 1e27 [nm^3]
        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_nm, 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::nm_to_m, 3));

        MJOLNIR_LOG_INFO("1 [m] = ", phys_type::m_to_length(), " [", length, "]");
        MJOLNIR_LOG_INFO("1 [", length, "] = ", phys_type::length_to_m(), " [m]");

        MJOLNIR_LOG_INFO("1 [L] = ", phys_type::L_to_volume(), " [", length, "^3]");
        MJOLNIR_LOG_INFO("1 [", length, "^3] = ", phys_type::volume_to_L(), " [L]");
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_units: unknown unit for length: `",
            toml::find(units, "length"), "here", {
            "expected value is one of the following.",
            "- \"nm\"",
            "- \"angstrom\""
            }));
    }

    MJOLNIR_LOG_INFO(u8"phys::kB = ", phys_type::kB(), " [", energy, "]");
    MJOLNIR_LOG_INFO(u8"phys::NA = ", phys_type::NA(), " [1/mol]");
    MJOLNIR_LOG_INFO(u8"phys::ε0 = ", phys_type::eps0(),
                     " [e^2 / (", energy, '*', length, ")]");

    return read_simulator<traitsT>(data);
}

} // mjolnir
#endif// MJOLNIR_READ_UNIT_SYSTEM_HPP
