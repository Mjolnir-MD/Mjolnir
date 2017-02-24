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

        std::string boundary("Unlimited");
        try{toml::get<toml::String>(iter->at("boundary"));}
        catch(std::out_of_range& except)
        {
            MJOLNIR_LOG_WARN("LocalForceField: boundary condition is not set.",
                    "UnlimitedBoundary is used");
        }
        MJOLNIR_LOG_INFO("boundary name", boundary);

        if(interaction == "BondLength")
        {
            for(auto para = params.cbegin(); para != params.cend(); ++para)
            {
                auto indices = toml::get<toml::Array<toml::Integer>>(
                        para->at("indices"));
                std::array<std::size_t, 2> idxs;
                idxs[0] = indices.at(0);
                idxs[1] = indices.at(1);

                lff.emplace_2body(read_2body_interaction<traitsT>(
                            interaction, boundary, potential, *para),
                        std::move(idxs));
            }
        }
        else if(interaction == "BondAngle")
        {
            for(auto para = params.cbegin(); para != params.cend(); ++para)
            {
                auto indices = toml::get<toml::Array<toml::Integer>>(
                        para->at("indices"));
                std::array<std::size_t, 3> idxs;
                idxs[0] = indices.at(0);
                idxs[1] = indices.at(1);
                idxs[2] = indices.at(2);

                lff.emplace_3body(read_3body_interaction<traitsT>(
                            interaction, boundary, potential, *para),
                        std::move(idxs));
            }
        }
        else if(interaction == "DihedralAngle")
        {
            for(auto para = params.cbegin(); para != params.cend(); ++para)
            {
                auto indices = toml::get<toml::Array<toml::Integer>>(
                        para->at("indices"));
                std::array<std::size_t, 4> idxs;
                idxs[0] = indices.at(0);
                idxs[1] = indices.at(1);
                idxs[2] = indices.at(2);
                idxs[3] = indices.at(3);
                lff.emplace_4body(read_4body_interaction<traitsT>(
                            interaction, boundary, potential, *para),
                        std::move(idxs));
            }
        }
        else
            throw std::runtime_error("unknown interaction: " + interaction);
    }
    MJOLNIR_LOG_DEBUG("read_local_force_field RETUENED");
    return lff;
}


} // mjolnir
#endif /*MJOLNIR_IO_TOML_READ_LOCAL_FORCEFIELD */
