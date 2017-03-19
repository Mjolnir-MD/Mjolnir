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
        const std::string interaction =
            toml::get<toml::String>(iter->at("interaction"));
        MJOLNIR_LOG_INFO("interaction name", interaction);

        const std::string potential =
            toml::get<toml::String>(iter->at("potential"));
        MJOLNIR_LOG_INFO("potential name", potential);

        std::string boundary("Unlimited");
        try{boundary = toml::get<toml::String>(iter->at("boundary"));}
        catch(std::out_of_range& except)
        {
            MJOLNIR_LOG_WARN("boundary setting not found.",
                             "UnlimitedBoundary is used.");
        }
        MJOLNIR_LOG_INFO("boundary name", boundary);

        const toml::Array<toml::Table> params =
            toml::get<toml::Array<toml::Table>>(iter->at("parameters"));
        MJOLNIR_LOG_INFO("parameter table size", params.size());

        lff.emplace(read_local_interaction<traitsT>(
                    interaction, potential, boundary, params));
    }
    MJOLNIR_LOG_DEBUG("read_local_force_field RETUENED");
    return lff;
}


} // mjolnir
#endif /*MJOLNIR_IO_TOML_READ_LOCAL_FORCEFIELD */
