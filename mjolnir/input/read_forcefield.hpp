#ifndef MJOLNIR_READ_FORCEFIELD
#define MJOLNIR_READ_FORCEFIELD
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/input/get_toml_value.hpp>
#include <mjolnir/input/read_interaction.hpp>

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
    const auto& ffs = detail::value_at(data, "forcefields", "<root>"
            ).cast<toml::value_t::Array>();

    if(ffs.size() <= N)
    {
        throw std::out_of_range("no enough forcefields: " + std::to_string(N));
    }
    const auto& ff = ffs.at(N).cast<toml::value_t::Table>();

    std::vector<toml::Table> fflocal;
    std::vector<toml::Table> ffglobal;

    if(ffs.count("local"))
    {
        fflocal  = toml::get<std::vector<toml::Table>>(ff.at("local"));
    }
    if(ffs.count("global"))
    {
        ffglobal = toml::get<std::vector<toml::Table>>(ff.at("global"));
    }
    return ForceField<traitsT>(
            read_local_forcefield<traitsT>(std::move(fflocal)),
            read_global_forcefield<traitsT>(std::move(ffglobal)));
}


} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD
