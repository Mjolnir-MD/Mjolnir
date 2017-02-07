#ifndef MJOLNIR_IO_TOML_READ_LOCAL_FORCEFIELD
#define MJOLNIR_IO_TOML_READ_LOCAL_FORCEFIELD
#include "read_local_potential.hpp"
#include "read_local_interaction.hpp"
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
LocalForceField<traitsT>
read_local_force_field(const toml::Array<toml::Table>& lffs)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_local_force_field CALLED");
    LocalForceField<traitsT> lff;
    for(auto iter = lffs.cbegin(); iter != lffs.cend(); ++iter)
    {
        const std::string potential =
            toml::get<toml::String>(iter->at("potential"));
        MJOLNIR_LOG_INFO("potential name", potential);

        const toml::Array<toml::Table> params =
            toml::get<toml::Array<toml::Table>>(iter->at("parameters"));
        MJOLNIR_LOG_INFO("parameter table size", params.size());

        const std::string interaction =
            toml::get<toml::String>(iter->at("interaction"));
        MJOLNIR_LOG_INFO("interaction name", interaction);

        const std::string boundary =
            toml::get<toml::String>(iter->at("boundary"));
        MJOLNIR_LOG_INFO("boundary name", boundary);

        if(interaction == "BondLength")
        {
            lff.emplace_2body(
                read_2body_interaction<traitsT>(interaction, boundary),
                read_local_potential_array<traitsT, 2>(potential, params));
        }
        else if(interaction == "BondAngle")
        {
            lff.emplace_3body(
                read_3body_interaction<traitsT>(interaction, boundary),
                read_local_potential_array<traitsT, 3>(potential, params));
        }
        else if(interaction == "DihedralAngle")
        {
            lff.emplace_4body(
                read_4body_interaction<traitsT>(interaction, boundary),
                read_local_potential_array<traitsT, 4>(potential, params));
        }
        else
            throw std::runtime_error("unknown interaction: " + interaction);
    }
    MJOLNIR_LOG_DEBUG("read_local_force_field RETUENED");
    return lff;
}


} // mjolnir
#endif /*MJOLNIR_IO_TOML_READ_LOCAL_FORCEFIELD */
