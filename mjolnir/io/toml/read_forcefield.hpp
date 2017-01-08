#ifndef MJOLNIR_TOML_READ_FORCEFIELD
#define MJOLNIR_TOML_READ_FORCEFIELD
#include <mjolnir/io/toml/read_global_forcefield.hpp>
#include <mjolnir/io/toml/read_local_forcefield.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
ForceField<traitsT>
read_force_field(const toml::Table& tab)
{
    LocalForceField<traitsT> local;
    try
    {
        const toml::Array<toml::Table> lff = 
            toml::get<toml::Array<toml::Table>>(tab.at("localforcefield"));
        LocalForceField<traitsT> tmp = read_local_force_field<traitsT>(lff);
        local = std::move(tmp);
    }
    catch(std::exception& except){}

    GlobalForceField<traitsT> global;
    try
    {
        const toml::Array<toml::Table> gff = 
            toml::get<toml::Array<toml::Table>>(tab.at("globalforcefield"));
        GlobalForceField<traitsT> tmp = read_global_force_field<traitsT>(gff);
        global = std::move(tmp);
    }
    catch(std::exception& except){}

    ForceField<traitsT> ff(std::move(local), std::move(global));
    return ff;
}

}//mjolnir
#endif/* MJOLNIR_TOML_READ_FORCEFIELD */
