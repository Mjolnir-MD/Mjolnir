#ifndef MJOLNIR_IO_TOML_READ_GLOBAL_POTENTIAL
#define MJOLNIR_IO_TOML_READ_GLOBAL_POTENTIAL
#include <mjolnir/potential/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/LennardJonesPotential.hpp>
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<GlobalPotentialBase<traitsT>>
read_excluded_volume(const toml::Table& potent)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_excluded_volume CALLED");

    toml::Array<toml::Table> param_table =
        toml::get<toml::Array<toml::Table>>(potent.at("parameters"));
    MJOLNIR_LOG_INFO("parameter table size", param_table.size());

    std::vector<typename ExcludedVolumePotential<traitsT>::parameter_type>
        params(param_table.size());
    for(auto iter = make_zip(param_table.cbegin(), params.begin());
            iter != make_zip(param_table.cend(), params.end()); ++iter)
    {
        const typename traitsT::real_type sigma =
            toml::get<toml::Float>(get<0>(iter)->at("sigma"));
        MJOLNIR_LOG_INFO("sigma", sigma);

        *get<1>(iter) = sigma;
    }
    const typename traitsT::real_type epsilon =
        toml::get<toml::Float>(potent.at("epsilon"));
    MJOLNIR_LOG_INFO("epsilon", epsilon);

    MJOLNIR_LOG_DEBUG("read_global_potential RETURNED");
    return make_unique<ExcludedVolumePotential<traitsT>>(
            epsilon, std::move(params));
}

template<typename traitsT>
std::unique_ptr<GlobalPotentialBase<traitsT>>
read_lennard_jones(const toml::Table& potent)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_lennard_jones CALLED");
    toml::Array<toml::Table> param_table =
        toml::get<toml::Array<toml::Table>>(potent.at("parameters"));
    MJOLNIR_LOG_INFO("parameter table size", param_table.size());

    std::vector<typename LennardJonesPotential<traitsT>::parameter_type>
        params(param_table.size());

    for(auto iter = make_zip(param_table.cbegin(), params.begin());
            iter != make_zip(param_table.cend(), params.end()); ++iter)
    {
        const typename traitsT::real_type sigma =
            toml::get<toml::Float>(get<0>(iter)->at("sigma"));
        const typename traitsT::real_type epsilon =
            toml::get<toml::Float>(get<0>(iter)->at("epsilon"));
        MJOLNIR_LOG_INFO("sigma", sigma, "epsilon", epsilon);

        *get<1>(iter) = std::make_pair(sigma, epsilon);
    }

    MJOLNIR_LOG_DEBUG("read_global_potential RETURNED");
    return make_unique<LennardJonesPotential<traitsT>>(
            std::move(params));
}

template<typename traitsT>
std::unique_ptr<GlobalPotentialBase<traitsT>>
read_global_potential(const std::string& name, const toml::Table& potent)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_global_potential CALLED");
    MJOLNIR_LOG_INFO("global potential name", name);

    if(name == "ExcludedVolume")
    {
        return read_excluded_volume<traitsT>(potent);
    }
    else if(name == "LennardJones")
    {
        return read_lennard_jones<traitsT>(potent);
    }
    else
        throw std::runtime_error("unknown potential: " + name);
}

} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_GLOBAL_POTENTIAL */
