#ifndef MJOLNIR_READ_UNIT_SYSTEM_HPP
#define MJOLNIR_READ_UNIT_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/core/constants.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_simulator.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<SimulatorBase> read_units(const toml::Table& data)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_units(const toml::Table& data), 0);
    using real_type = typename traitsT::real_type;
    using phys_type = physics::constants<real_type>;
    using math_type = math::constants<real_type>;
    using unit_type = unit::constants<real_type>;

    const auto& units = get_toml_value<toml::Table>(data, "units", "<root>");

    const auto& energy = get_toml_value<std::string>(units, "energy", "[units]");
    const auto& length = get_toml_value<std::string>(units, "length", "[units]");

    MJOLNIR_LOG_INFO("unit of energy : [", energy, ']');
    MJOLNIR_LOG_INFO("unit of length : [", length, ']');

    if(energy == "kcal/mol")
    {
        // kB [J/K] -> [kcal/mol/K] by * (J to cal) * 1e-3 * (/mol)
        phys_type::kB   *= (unit_type::J_to_cal / 1000.0) *
                           unit_type::avogadro_constant;
        // eps0 [F/m] == [C^2/J/m] -> [C^2/(kcal/mol)/m]
        phys_type::eps0 *= (1000.0 / unit_type::J_to_cal) / 
                           unit_type::avogadro_constant;
    }
    else if(energy == "kJ/mol")
    {
        // kB [J/K] -> [kJ/mol/K]
        phys_type::kB   *= 1e-3 * unit_type::avogadro_constant;
        // eps0 [F/m] == [C^2/J/m] -> [C^2/kJ/mol/m]
        phys_type::eps0 *= 1e+3 / unit_type::avogadro_constant;
    }
    else
    {
        throw_exception<std::runtime_error>("mjolnir::read_units: unknown unit "
            "for energy: `", energy, "`. `kcal/mol`, `kJ/mol` are allowed");
    }

    // until here, SI `m` are used as length unit.

    if(length == "angstrom" || length == "Å")
    {
        // eps0 [C^2/Energy/m] -> [C^2/Energy/Angstrom]
        phys_type::eps0 *= (1.0 / unit_type::m_to_angstrom);
    }
    else if(length == "nm")
    {
        // eps0 [C^2/Energy/m] -> [C^2/Energy/nm]
        phys_type::eps0 *= (1.0 / unit_type::m_to_nm);
    }
    else
    {
        throw_exception<std::runtime_error>("mjolnir::read_units: unknown unit "
            "for length: `", length, "`. `angstrom`, `nm` are allowed");
    }

    MJOLNIR_LOG_INFO(u8"phys::kB = ", phys_type::kB,   " [", energy, "]");
    MJOLNIR_LOG_INFO(u8"phys::NA = ", phys_type::NA,   " [1/mol]");
    MJOLNIR_LOG_INFO(u8"phys::e  = ", phys_type::e,    " [C]");
    MJOLNIR_LOG_INFO(u8"phys::ε0 = ", phys_type::eps0,
                       " [C^2 / (", energy, '*', length, ")]");

    return read_simulator<traitsT>(data);
}

} // mjolnir
#endif// MJOLNIR_READ_UNIT_SYSTEM_HPP
