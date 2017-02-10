#ifndef MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD
#define MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD
#include "read_global_potential.hpp"
#include "read_global_interaction.hpp"
#include <mjolnir/core/GlobalForceField.hpp>
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
GlobalForceField<traitsT>
read_global_force_field(const toml::Array<toml::Table>& gffs)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_global_force_field CALLED");
    GlobalForceField<traitsT> gff;

    for(auto iter = gffs.cbegin(); iter != gffs.cend(); ++iter)
    {
        const std::string potential =
            toml::get<toml::String>(iter->at("potential"));
        MJOLNIR_LOG_INFO("potential name read", potential);

        std::unique_ptr<GlobalPotentialBase<traitsT>> pot =
            read_global_potential<traitsT>(potential, *iter);

        const std::string interaction =
            toml::get<toml::String>(iter->at("interaction"));
        MJOLNIR_LOG_INFO("interaction name read", interaction);

        std::unique_ptr<GlobalInteractionBase<traitsT>> inter =
            read_global_interaction<traitsT>(interaction, pot, *iter);

        gff.emplace(std::move(inter), std::move(pot));
    }
    MJOLNIR_LOG_DEBUG("read_global_force_field RETURNED");
    return gff;
}

} // mjolnir
#endif/* MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD */
