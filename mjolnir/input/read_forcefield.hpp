#ifndef MJOLNIR_READ_FORCEFIELD
#define MJOLNIR_READ_FORCEFIELD
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include "read_interaction.hpp"

namespace mjolnir
{

template<typename traitsT>
LocalForceField<traitsT>
read_local_forcefield(std::vector<toml::Table> interactions)
{
    LocalForceField<traitsT> lff;
    //XXX: TODO?
    // by using typeid(), merge interactions if the types are same.
    for(const auto& interaction : interactions)
    {
        lff.emplace(read_local_interaction<traitsT>(interaction));
    }
    return lff;
}


template<typename traitsT>
GlobalForceField<traitsT>
read_global_forcefield(std::vector<toml::Table> interactions)
{
    GlobalForceField<traitsT> gff;
    //XXX: TODO?
    // by using typeid(), merge interactions if the types are same.
    for(const auto& interaction : interactions)
    {
        gff.emplace(read_global_interaction<traitsT>(interaction));
    }
    return gff;
}

template<typename traitsT>
ForceField<traitsT>
read_forcefield(const toml::Table& data, std::size_t N)
{
    const auto& ffs = data.at("forcefields").cast<toml::value_t::Array>();
    if(ffs.size() <= N)
        throw std::out_of_range("no enough forcefields: " + std::to_string(N));
    const auto& ff = ffs.at(N).cast<toml::value_t::Table>();

    return ForceField<traitsT>(
        read_local_forcefield<traitsT>(
            toml::get<std::vector<toml::Table>>(ff.at("local"))),
        read_global_forcefield<traitsT>(
            toml::get<std::vector<toml::Table>>(ff.at("global"))));
}


} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD
